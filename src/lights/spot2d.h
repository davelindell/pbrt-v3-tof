
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LIGHTS_SPOT2D_H
#define PBRT_LIGHTS_SPOT2D_H

// lights/spot.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"

namespace pbrt {

// SpotLight2d Declarations
class SpotLight2d : public Light {
  public:
    // SpotLight2d Public Methods
    SpotLight2d(const Transform &l2w, const Point3f &from,
              const Transform &LightToWorld, const MediumInterface &m,
              const Spectrum &I, Float totalWidth, Float falloffStart, int confocal);
    SpotLight2d(SpotLight2d &light);
    Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi,
                       Float *pdf, VisibilityTester *vis) const;
    Float Falloff(const Vector3f &w) const;
    Spectrum Power() const;
    Float Pdf_Li(const Interaction &, const Vector3f &) const;
    Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                       Ray *ray, Normal3f *nLight, Float *pdfPos,
                       Float *pdfDir) const;
    void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
                Float *pdfDir) const;
    void AdjustDirection(Vector3f &dir);
    std::shared_ptr<Light> Clone();
    std::shared_ptr<SpotLight2d> doClone();

  private:
    // SpotLight2d Private Data
    const Point3f pLight;
    const Spectrum I;
    const Float cosTotalWidth, cosFalloffStart;
    Transform LToW;
    Transform WToL;
    const Transform l2w;
    const Point3f from;
};

std::shared_ptr<SpotLight2d> CreateSpotLight2d(const Transform &l2w,
                                           const Medium *medium,
                                           const ParamSet &paramSet);

}  // namespace pbrt

#endif  // PBRT_LIGHTS_SPOT_H
