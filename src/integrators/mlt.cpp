
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/mlt.cpp*
#include "integrators/mlt.h"
#include "integrators/bdpt.h"
#include "scene.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "camera.h"
#include "stats.h"
#include "filters/box.h"
#include "paramset.h"
#include "sampling.h"
#include "progressreporter.h"

namespace pbrt {

STAT_PERCENT("Integrator/Acceptance rate", acceptedMutations, totalMutations);

// MLTSampler Constants
static const int cameraStreamIndex = 0;
static const int lightStreamIndex = 1;
static const int connectionStreamIndex = 2;
static const int nSampleStreams = 3;

// MLTSampler Method Definitions
Float MLTSampler::Get1D() {
    ProfilePhase _(Prof::GetSample);
    int index = GetNextIndex();
    EnsureReady(index);
    return X[index].value;
}

Point2f MLTSampler::Get2D() { return {Get1D(), Get1D()}; }

std::unique_ptr<Sampler> MLTSampler::Clone(int seed) {
    LOG(FATAL) << "MLTSampler::Clone() is not implemented";
    return nullptr;
}

void MLTSampler::StartIteration() {
    currentIteration++;
    largeStep = rng.UniformFloat() < largeStepProbability;
}

void MLTSampler::Accept() {
    if (largeStep) lastLargeStepIteration = currentIteration;
}

void MLTSampler::EnsureReady(int index) {
    // Enlarge _MLTSampler::X_ if necessary and get current $\VEC{X}_i$
    if (index >= X.size()) X.resize(index + 1);
    PrimarySample &Xi = X[index];

    // Reset $\VEC{X}_i$ if a large step took place in the meantime
    if (Xi.lastModificationIteration < lastLargeStepIteration) {
        Xi.value = rng.UniformFloat();
        Xi.lastModificationIteration = lastLargeStepIteration;
    }

    // Apply remaining sequence of mutations to _sample_
    Xi.Backup();
    if (largeStep) {
        Xi.value = rng.UniformFloat();
    } else {
        int64_t nSmall = currentIteration - Xi.lastModificationIteration;
        // Apply _nSmall_ small step mutations

        // Sample the standard normal distribution $N(0, 1)$
        Float normalSample = Sqrt2 * ErfInv(2 * rng.UniformFloat() - 1);

        // Compute the effective standard deviation and apply perturbation to
        // $\VEC{X}_i$
        Float effSigma = sigma * std::sqrt((Float)nSmall);
        Xi.value += normalSample * effSigma;
        Xi.value -= std::floor(Xi.value);
    }
    Xi.lastModificationIteration = currentIteration;
}

void MLTSampler::Reject() {
    for (auto &Xi : X)
        if (Xi.lastModificationIteration == currentIteration) Xi.Restore();
    --currentIteration;
}

void MLTSampler::StartStream(int index) {
    CHECK_LT(index, streamCount);
    streamIndex = index;
    sampleIndex = 0;
}

// MLT Method Definitions
Spectrum MLTIntegrator::L(const Scene &scene, MemoryArena &arena,
                          const std::unique_ptr<Distribution1D> &lightDistr,
                          const std::unordered_map<const Light *, size_t> &lightToIndex,
                          MLTSampler &sampler, int depth, Point2f *pRaster) {
    sampler.StartStream(cameraStreamIndex);
    // Determine the number of available strategies and pick a specific one
    int s, t, nStrategies;
    if (depth == 0) {
        nStrategies = 1;
        s = 0;
        t = 2;
    } else {
        nStrategies = depth + 2;
        s = std::min((int)(sampler.Get1D() * nStrategies), nStrategies - 1);
        t = nStrategies - s;
    }

    // Generate a camera subpath with exactly _t_ vertices
    Vertex *cameraVertices = arena.Alloc<Vertex>(t);
    Bounds2f sampleBounds = (Bounds2f)camera->film->GetSampleBounds();
    *pRaster = sampleBounds.Lerp(sampler.Get2D());
    if (GenerateCameraSubpath(scene, sampler, arena, t, *camera, *pRaster,
                              cameraVertices) != t)
        return Spectrum(0.f);

    // Generate a light subpath with exactly _s_ vertices
    sampler.StartStream(lightStreamIndex);
    Vertex *lightVertices = arena.Alloc<Vertex>(s);
    if (GenerateLightSubpath(scene, sampler, arena, s, cameraVertices[0].time(),
                             *lightDistr, lightToIndex, lightVertices) != s)
        return Spectrum(0.f);

    // Execute connection strategy and return the radiance estimate
    sampler.StartStream(connectionStreamIndex);
    return ConnectBDPT(scene, lightVertices, cameraVertices, s, t, *lightDistr,
                       lightToIndex, *camera, sampler, pRaster) *
           nStrategies;
}

void MLTIntegrator::Render(const Scene &scene) {
    std::unique_ptr<Distribution1D> lightDistr =
        ComputeLightPowerDistribution(scene);

    // Compute a reverse mapping from light pointers to offsets into the
    // scene lights vector (and, equivalently, offsets into
    // lightDistr). Added after book text was finalized; this is critical
    // to reasonable performance with 100s+ of light sources.
    std::unordered_map<const Light *, size_t> lightToIndex;
    for (size_t i = 0; i < scene.lights.size(); ++i)
        lightToIndex[scene.lights[i].get()] = i;

    // Generate bootstrap samples and compute normalization constant $b$
    int nBootstrapSamples = nBootstrap * (maxDepth + 1);
    std::vector<Float> bootstrapWeights(nBootstrapSamples, 0);
    if (scene.lights.size() > 0) {
        ProgressReporter progress(nBootstrap / 256,
                                  "Generating bootstrap paths");
        std::vector<MemoryArena> bootstrapThreadArenas(MaxThreadIndex());
        int chunkSize = Clamp(nBootstrap / 128, 1, 8192);
        ParallelFor([&](int i) {
            // Generate _i_th bootstrap sample
            MemoryArena &arena = bootstrapThreadArenas[ThreadIndex];
            for (int depth = 0; depth <= maxDepth; ++depth) {
                int rngIndex = i * (maxDepth + 1) + depth;
                MLTSampler sampler(mutationsPerPixel, rngIndex, sigma,
                                   largeStepProbability, nSampleStreams);
                Point2f pRaster;
                bootstrapWeights[rngIndex] =
                    L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pRaster).y();
                arena.Reset();
            }
            if ((i + 1) % 256 == 0) progress.Update();
        }, nBootstrap, chunkSize);
        progress.Done();
    }
    Distribution1D bootstrap(&bootstrapWeights[0], nBootstrapSamples);
    Float b = bootstrap.funcInt * (maxDepth + 1);

    // Run _nChains_ Markov chains in parallel
    Film &film = *camera->film;
    int64_t nTotalMutations =
        (int64_t)mutationsPerPixel * (int64_t)film.GetSampleBounds().Area();
    if (scene.lights.size() > 0) {
        const int progressFrequency = 32768;
        ProgressReporter progress(nTotalMutations / progressFrequency,
                                  "Rendering");
        ParallelFor([&](int i) {
            int64_t nChainMutations =
                std::min((i + 1) * nTotalMutations / nChains, nTotalMutations) -
                i * nTotalMutations / nChains;
            // Follow {i}th Markov chain for _nChainMutations_
            MemoryArena arena;

            // Select initial state from the set of bootstrap samples
            RNG rng(i);
            int bootstrapIndex = bootstrap.SampleDiscrete(rng.UniformFloat());
            int depth = bootstrapIndex % (maxDepth + 1);

            // Initialize local variables for selected state
            MLTSampler sampler(mutationsPerPixel, bootstrapIndex, sigma,
                               largeStepProbability, nSampleStreams);
            Point2f pCurrent;
            Spectrum LCurrent =
                L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pCurrent);

            // Run the Markov chain for _nChainMutations_ steps
            for (int64_t j = 0; j < nChainMutations; ++j) {
                sampler.StartIteration();
                Point2f pProposed;
                Spectrum LProposed =
                    L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pProposed);
                // Compute acceptance probability for proposed sample
                Float accept = std::min((Float)1, LProposed.y() / LCurrent.y());

                // Splat both current and proposed samples to _film_
                if (accept > 0)
                    film.AddSplat(pProposed,
                                  LProposed * accept / LProposed.y());
                film.AddSplat(pCurrent, LCurrent * (1 - accept) / LCurrent.y());

                // Accept or reject the proposal
                if (rng.UniformFloat() < accept) {
                    pCurrent = pProposed;
                    LCurrent = LProposed;
                    sampler.Accept();
                    ++acceptedMutations;
                } else
                    sampler.Reject();
                ++totalMutations;
                if ((i * nTotalMutations / nChains + j) % progressFrequency ==
                    0)
                    progress.Update();
                arena.Reset();
            }
        }, nChains);
        progress.Done();
    }

    // Store final image computed with MLT
    camera->film->WriteImage(b / mutationsPerPixel);
}

//int GenerateCameraSubpath(const Scene &scene, Sampler &sampler,
//                          MemoryArena &arena, int maxDepth,
//                          const Camera &camera, const Point2f &pFilm,
//                          Vertex *path) {
//    if (maxDepth == 0) return 0;
//    ProfilePhase _(Prof::BDPTGenerateSubpath);
//    // Sample initial ray for camera subpath
//    CameraSample cameraSample;
//    cameraSample.pFilm = pFilm;
//    cameraSample.time = sampler.Get1D();
//    cameraSample.pLens = sampler.Get2D();
//    RayDifferential ray;
//    Spectrum beta = camera.GenerateRayDifferential(cameraSample, &ray);
//    ray.ScaleDifferentials(1 / std::sqrt(sampler.samplesPerPixel));
//
//    // adjust spotlight to point along the camera ray direction
//    for (size_t i = 0; i < scene.lights.size(); ++i)
//        scene.lights[i].get()->AdjustDirection(ray.d);
//
//    // Generate first vertex on camera subpath and start random walk
//    Float pdfPos, pdfDir;
//    path[0] = Vertex::CreateCamera(&camera, ray, beta);
//    camera.Pdf_We(ray, &pdfPos, &pdfDir);
//    VLOG(2) << "Starting camera subpath. Ray: " << ray << ", beta " << beta
//            << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;
//    return RandomWalk(scene, ray, sampler, arena, beta, pdfDir, maxDepth - 1,
//                      TransportMode::Radiance, path + 1) +
//           1;
//}
//
//int GenerateLightSubpath(
//    const Scene &scene, Sampler &sampler, MemoryArena &arena, int maxDepth,
//    Float time, const Distribution1D &lightDistr,
//    const std::unordered_map<const Light *, size_t> &lightToIndex,
//    Vertex *path) {
//    if (maxDepth == 0) return 0;
//    ProfilePhase _(Prof::BDPTGenerateSubpath);
//    // Sample initial ray for light subpath
//    Float lightPdf;
//    int lightNum = lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
//    const std::shared_ptr<Light> &light = scene.lights[lightNum];
//    RayDifferential ray;
//    Normal3f nLight;
//    Float pdfPos, pdfDir;
//    Spectrum Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(), time, &ray,
//                                   &nLight, &pdfPos, &pdfDir);
//    if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) return 0;
//
//    // Generate first vertex on light subpath and start random walk
//    path[0] =
//        Vertex::CreateLight(light.get(), ray, nLight, Le, pdfPos * lightPdf);
//    Spectrum beta = Le * AbsDot(nLight, ray.d) / (lightPdf * pdfPos * pdfDir);
//    VLOG(2) << "Starting light subpath. Ray: " << ray << ", Le " << Le <<
//        ", beta " << beta << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;
//    int nVertices =
//        RandomWalk(scene, ray, sampler, arena, beta, pdfDir, maxDepth - 1,
//                   TransportMode::Importance, path + 1);
//
//    // Correct subpath sampling densities for infinite area lights
//    if (path[0].IsInfiniteLight()) {
//        // Set spatial density of _path[1]_ for infinite area light
//        if (nVertices > 0) {
//            path[1].pdfFwd = pdfPos;
//            if (path[1].IsOnSurface())
//                path[1].pdfFwd *= AbsDot(ray.d, path[1].ng());
//        }
//
//        // Set spatial density of _path[0]_ for infinite area light
//        path[0].pdfFwd =
//            InfiniteLightDensity(scene, lightDistr, lightToIndex, ray.d);
//    }
//    return nVertices + 1;
//}
//
//Spectrum ConnectBDPT(
//    const Scene &scene, Vertex *lightVertices, Vertex *cameraVertices, int s,
//    int t, const Distribution1D &lightDistr,
//    const std::unordered_map<const Light *, size_t> &lightToIndex,
//    const Camera &camera, Sampler &sampler, Point2f *pRaster,
//    Float *misWeightPtr) {
//    ProfilePhase _(Prof::BDPTConnectSubpaths);
//    Spectrum L(0.f);
//    // Ignore invalid connections related to infinite area lights
//    if (t > 1 && s != 0 && cameraVertices[t - 1].type == VertexType::Light)
//        return Spectrum(0.f);
//
//    // Perform connection and write contribution to _L_
//    Vertex sampled;
//    if (s == 0) {
//        // Interpret the camera subpath as a complete path
//        const Vertex &pt = cameraVertices[t - 1];
//        if (pt.IsLight()) L = pt.Le(scene, cameraVertices[t - 2]) * pt.beta;
//        DCHECK(!L.HasNaNs());
//        float distance = 0;
//        for (int i = 0; i < t-1; ++i) {
//            distance += Distance(cameraVertices[i].p(), cameraVertices[i+1].p()); 
//        }
//        L = L.AddTimeOfFlight(distance);
//
//    } else if (t == 1) {
//        // Sample a point on the camera and connect it to the light subpath
//        const Vertex &qs = lightVertices[s - 1];
//        if (qs.IsConnectible()) {
//            VisibilityTester vis;
//            Vector3f wi;
//            Float pdf;
//            Spectrum Wi = camera.Sample_Wi(qs.GetInteraction(), sampler.Get2D(),
//                                           &wi, &pdf, pRaster, &vis);
//            if (pdf > 0 && !Wi.IsBlack()) {
//                // Initialize dynamically sampled vertex and _L_ for $t=1$ case
//                sampled = Vertex::CreateCamera(&camera, vis.P1(), Wi / pdf);
//                L = qs.beta * qs.f(sampled, TransportMode::Importance) * sampled.beta;
//                if (qs.IsOnSurface()) L *= AbsDot(wi, qs.ns());
//                DCHECK(!L.HasNaNs());
//                // Only check visibility after we know that the path would
//                // make a non-zero contribution.
//                if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
//
//                // calculate distance camera -> light
//                float distance = 0;
//                for (int i = 0; i < s-1; ++i) {
//                    distance += Distance(lightVertices[i].p(), lightVertices[i+1].p()); 
//                }
//                distance += Distance(lightVertices[s-1].p(), sampled.p());
//                L = L.AddTimeOfFlight(distance);
//            }
//        }
//    } else if (s == 1) {
//        // Sample a point on a light and connect it to the camera subpath
//        const Vertex &pt = cameraVertices[t - 1];
//        if (pt.IsConnectible()) {
//            Float lightPdf;
//            VisibilityTester vis;
//            Vector3f wi;
//            Float pdf;
//            int lightNum =
//                lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
//            const std::shared_ptr<Light> &light = scene.lights[lightNum];
//            Spectrum lightWeight = light->Sample_Li(
//                pt.GetInteraction(), sampler.Get2D(), &wi, &pdf, &vis);
//            if (pdf > 0 && !lightWeight.IsBlack()) {
//                EndpointInteraction ei(vis.P1(), light.get());
//                sampled =
//                    Vertex::CreateLight(ei, lightWeight / (pdf * lightPdf), 0);
//                sampled.pdfFwd =
//                    sampled.PdfLightOrigin(scene, pt, lightDistr, lightToIndex);
//                L = pt.beta * pt.f(sampled, TransportMode::Radiance) * sampled.beta;
//                if (pt.IsOnSurface()) L *= AbsDot(wi, pt.ns());
//                // Only check visibility if the path would carry radiance.
//                if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
//
//                // calculate distance camera -> light
//                float distance = 0;
//                for (int i = 0; i < t-1; ++i) {
//                    distance += Distance(cameraVertices[i].p(), cameraVertices[i+1].p()); 
//                }
//                distance += Distance(cameraVertices[t-1].p(), sampled.p());
//                L = L.AddTimeOfFlight(distance);
//
//                //L = L.AddTimeOfFlight(Distance(sampled.p(), pt.p()));
//            }
//        }
//    } else {
//        // Handle all other bidirectional connection cases
//        const Vertex &qs = lightVertices[s - 1], &pt = cameraVertices[t - 1];
//        if (qs.IsConnectible() && pt.IsConnectible()) {
//            L = qs.beta * qs.f(pt, TransportMode::Importance) * pt.f(qs, TransportMode::Radiance) * pt.beta;
//            VLOG(2) << "General connect s: " << s << ", t: " << t <<
//                " qs: " << qs << ", pt: " << pt << ", qs.f(pt): " << qs.f(pt, TransportMode::Importance) <<
//                ", pt.f(qs): " << pt.f(qs, TransportMode::Radiance) << ", G: " << G(scene, sampler, qs, pt) <<
//                ", dist^2: " << DistanceSquared(qs.p(), pt.p());
//            if (!L.IsBlack()) L *= G(scene, sampler, qs, pt);
//
//            float distance = 0;
//            for (int i = 0; i < t-1; ++i) {
//                distance += Distance(cameraVertices[i].p(), cameraVertices[i+1].p()); 
//            }
//            distance += Distance(cameraVertices[t-1].p(), lightVertices[s-1].p());
//            for (int i = 0; i < s-1; ++i) {
//                distance += Distance(lightVertices[i].p(), lightVertices[i+1].p()); 
//            }
//            L = L.AddTimeOfFlight(distance);
//        }
//    }
//
//    ++totalPaths;
//    if (L.IsBlack()) ++zeroRadiancePaths;
//    ReportValue(pathLength, s + t - 2);
//
//    // Compute MIS weight for connection strategy
//    Float misWeight =
//        L.IsBlack() ? 0.f : MISWeight(scene, lightVertices, cameraVertices,
//                                      sampled, s, t, lightDistr, lightToIndex);
//    VLOG(2) << "MIS weight for (s,t) = (" << s << ", " << t << ") connection: "
//            << misWeight;
//    DCHECK(!std::isnan(misWeight));
//    L *= misWeight;
//    if (misWeightPtr) *misWeightPtr = misWeight;
//    return L;
//}

MLTIntegrator *CreateMLTIntegrator(const ParamSet &params,
                                   std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
    int64_t nChains = params.FindOneInt("chains", 1000);
    int mutationsPerPixel = params.FindOneInt("mutationsperpixel", 100);
    Float largeStepProbability =
        params.FindOneFloat("largestepprobability", 0.3f);
    Float sigma = params.FindOneFloat("sigma", .01f);
    if (PbrtOptions.quickRender) {
        mutationsPerPixel = std::max(1, mutationsPerPixel / 16);
        nBootstrap = std::max(1, nBootstrap / 16);
    }
    return new MLTIntegrator(camera, maxDepth, nBootstrap, nChains,
                             mutationsPerPixel, sigma, largeStepProbability);
}

}  // namespace pbrt
