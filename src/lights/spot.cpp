
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


// lights/spot.cpp*
#include "lights/spot.h"
#include "paramset.h"
#include "sampling.h"
#include "reflection.h"
#include "stats.h"

namespace pbrt {

// SpotLight Method Definitions
SpotLight::SpotLight(const Transform &l2w,
                     const Point3f &from, 
                     const Transform &LightToWorld,
                     const MediumInterface &mediumInterface, const Spectrum &I,
                     Float totalWidth, Float falloffStart, int confocal)
    : Light((int)LightFlags::DeltaPosition, LightToWorld, mediumInterface, 1, confocal),
      l2w(l2w),
      LToW(LightToWorld),
      WToL(Inverse(LightToWorld)),
      from(from),
      pLight(LightToWorld(Point3f(0, 0, 0))),
      I(I),
      cosTotalWidth(std::cos(Radians(totalWidth))),
      cosFalloffStart(std::cos(Radians(falloffStart))) {}

SpotLight::SpotLight(SpotLight &light)
    : Light(light),
      l2w(light.l2w),
      LToW(light.LightToWorld),
      WToL(Inverse(light.LightToWorld)),
      from(light.from),
      pLight(light.LightToWorld(Point3f(0, 0, 0))),
      I(light.I),
      cosTotalWidth(light.cosTotalWidth),
      cosFalloffStart(light.cosFalloffStart) {}

Spectrum SpotLight::Sample_Li(const Interaction &ref, const Point2f &u,
                              Vector3f *wi, Float *pdf,
                              VisibilityTester *vis) const {
    ProfilePhase _(Prof::LightSample);
    *wi = Normalize(pLight - ref.p);
    *pdf = 1.f;
    *vis =
        VisibilityTester(ref, Interaction(pLight, ref.time, mediumInterface));
    return I * Falloff(-*wi) / DistanceSquared(pLight, ref.p);
}

Float SpotLight::Falloff(const Vector3f &w) const {
    Vector3f wl = Normalize(WToL(w));
    Float cosTheta = wl.z;
    if (cosTheta < cosTotalWidth) return 0;
    if (cosTheta >= cosFalloffStart) return 1;
    // Compute falloff inside spotlight cone
    Float delta =
        (cosTheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);
    return (delta * delta) * (delta * delta);
}

Spectrum SpotLight::Power() const {
    return I * 2 * Pi * (1 - .5f * (cosFalloffStart + cosTotalWidth));
}

Float SpotLight::Pdf_Li(const Interaction &, const Vector3f &) const {
    return 0.f;
}

Spectrum SpotLight::Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                              Ray *ray, Normal3f *nLight, Float *pdfPos,
                              Float *pdfDir) const {
    ProfilePhase _(Prof::LightSample);
    Vector3f w = UniformSampleCone(u1, cosTotalWidth);
    *ray = Ray(pLight, LToW(w), Infinity, time, mediumInterface.inside);
    *nLight = (Normal3f)ray->d;
    *pdfPos = 1;
    *pdfDir = UniformConePdf(cosTotalWidth);
    return I * Falloff(ray->d);
}

void SpotLight::Pdf_Le(const Ray &ray, const Normal3f &, Float *pdfPos,
                       Float *pdfDir) const {
    ProfilePhase _(Prof::LightPdf);
    *pdfPos = 0;
    *pdfDir = (CosTheta(WToL(ray.d)) >= cosTotalWidth)
                  ? UniformConePdf(cosTotalWidth)
                  : 0;
}

void SpotLight::AdjustDirection(Vector3f &dir)
{
Vector3f du, dv;
CoordinateSystem(dir, &du, &dv);
Transform dirToZ =
    Transform(Matrix4x4(du.x, du.y, du.z, 0., dv.x, dv.y, dv.z, 0., dir.x,
                        dir.y, dir.z, 0., 0, 0, 0, 1.));
LToW =
    l2w * Translate(Vector3f(from.x, from.y, from.z)) * Inverse(dirToZ);
WToL = Inverse(LToW); 
}

std::shared_ptr<SpotLight> SpotLight::doClone() {
    return std::make_shared<SpotLight>(*this);
}

std::shared_ptr<Light> SpotLight::Clone() {
    return doClone();
}

std::shared_ptr<SpotLight> CreateSpotLight(const Transform &l2w,
                                           const Medium *medium,
                                           const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Float coneangle = paramSet.FindOneFloat("coneangle", 30.);
    Float conedelta = paramSet.FindOneFloat("conedeltaangle", 5.);
    int confocal = paramSet.FindOneInt("confocal", 0);
    // Compute spotlight world to light transformation
    Point3f from = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
    Point3f to = paramSet.FindOnePoint3f("to", Point3f(0, 0, 1));
    Vector3f dir = Normalize(to - from);
    Vector3f du, dv;
    CoordinateSystem(dir, &du, &dv);
    Transform dirToZ =
        Transform(Matrix4x4(du.x, du.y, du.z, 0., dv.x, dv.y, dv.z, 0., dir.x,
                            dir.y, dir.z, 0., 0, 0, 0, 1.));
    Transform light2world =
        l2w * Translate(Vector3f(from.x, from.y, from.z)) * Inverse(dirToZ);
    return std::make_shared<SpotLight>(l2w, from, light2world, medium, I * sc, coneangle,
                                       coneangle - conedelta, confocal);
}

}  // namespace pbrt
