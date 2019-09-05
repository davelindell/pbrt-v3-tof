
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


// cameras/orthographic.cpp*
#include "cameras/orthographic2d.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "stats.h"
#include "light.h"

namespace pbrt {

// OrthographicCamera Definitions
Float OrthographicCamera2d::GenerateRay(const CameraSample &sample,
                                      Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, 0, 0);
    Point3f pCamera = RasterToCamera(pFilm);
    *ray = Ray(pCamera, Vector3f(0, 0, 1));
    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
}

Float OrthographicCamera2d::GenerateRayDifferential(const CameraSample &sample,
                                                  RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute main orthographic viewing ray

    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, 0, 0);
    Point3f pCamera = RasterToCamera(pFilm);
    *ray = RayDifferential(pCamera, Vector3f(0, 0, 1));

    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    // Compute ray differentials for _OrthographicCamera_
    if (lensRadius > 0) {
        // Compute _OrthographicCamera_ ray differentials accounting for lens

        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Float ft = focalDistance / ray->d.z;

        Point3f pFocus = pCamera + dxCamera + (ft * Vector3f(0, 0, 1));
        ray->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

        pFocus = pCamera + dyCamera + (ft * Vector3f(0, 0, 1));
        ray->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
    } else {
        ray->rxOrigin = ray->o + dxCamera;
        ray->ryOrigin = ray->o + dyCamera;
        ray->rxDirection = ray->ryDirection = ray->d;
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->hasDifferentials = true;
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
}

OrthographicCamera2d *CreateOrthographicCamera2d(const ParamSet &params,
                                             const AnimatedTransform &cam2world,
                                             Film *film, const Medium *medium) {
    // convert to film2d
    Film2d *film2d = new Film2d(film);
    delete film;

    // Extract common camera parameters from _ParamSet_
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6f);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film2d->fullResolution.x) / Float(film2d->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    return new OrthographicCamera2d(cam2world, screen, shutteropen, shutterclose,
                                  lensradius, focaldistance, film2d, medium);
}

Spectrum OrthographicCamera2d::We(const Ray &ray, Point2f *pRaster2) const {
    // Interpolate camera matrix and check if $\w{}$ is forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) return 0;

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return raster position if requested
    if (pRaster2) *pRaster2 = Point2f(pRaster.x, pRaster.y);

    // Return zero importance for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y)
        return 0;

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;

    // Return importance for point on image plane
    Float cos2Theta = cosTheta * cosTheta;
    //return Spectrum(1 / (lensArea * cos2Theta * cosTheta));
    return Spectrum(1 / (4*Pi));
}

void OrthographicCamera2d::Pdf_We(const Ray &ray, Float *pdfPos,
                               Float *pdfDir) const {
    // Interpolate camera matrix and fail if $\w{}$ is not forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return zero probability for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdfPos = 1 / lensArea;
    //*pdfDir = 1 / (cosTheta * cosTheta);
    *pdfDir = 1 / (4 * Pi);
}

Spectrum OrthographicCamera2d::Sample_Wi(const Interaction &ref, const Point2f &u,
                                      Vector3f *wi, Float *pdf,
                                      Point2f *pRaster,
                                      VisibilityTester *vis) const {
    // Uniformly sample a lens interaction _lensIntr_
    //Point2f pLens = lensRadius * ConcentricSampleDisk(u);
    Point2f pLens = ConcentricSampleDisk(u);
    Point3f pLensWorld = CameraToWorld(ref.time, Point3f(pLens.x, 0, 0));
    Interaction lensIntr(pLensWorld, ref.time, medium);
    lensIntr.n = Normal3f(CameraToWorld(ref.time, Vector3f(0, 0, 1)));

    // Populate arguments and compute the importance value
    *vis = VisibilityTester(ref, lensIntr);
    *wi = lensIntr.p - ref.p;
    Float dist = wi->Length();
    *wi /= dist;

    // Compute PDF for importance arriving at _ref_

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdf = (dist) / (AbsDot(lensIntr.n, *wi) * lensArea);
    //return We(lensIntr.SpawnRay(-*wi), pRaster);
    return Spectrum(1 / (4*Pi));
}


}  // namespace pbrt
