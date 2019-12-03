include_directories(.)
set(CRKSPH_inst
    CRKSPHHydroBase
    ContinuityVolumePolicy
    SolidCRKSPHHydroBase
    computeCRKSPHCorrections
    computeCRKSPHEvaluation
    computeCRKSPHIntegral
    computeCRKSPHMoments
    computeCRKSPHSumMassDensity
    computeCRKSPHSumVolume
    computeHVolumes
    computeHullVolumes
    computeSolidCRKSPHSumMassDensity
    detectSurface
    #    flagSurfaceNeighbors
    gradientCRKSPH
    interpolateCRKSPH
    zerothOrderSurfaceCorrections
   )


set(CRKSPH_sources
    centerOfMass.cc
    computeVoronoiVolume1d.cc
    computeVoronoiVolume2d.cc
    computeVoronoiVolume3d.cc
   )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHHydroBaseRZ.cc
       SolidCRKSPHHydroBaseRZ.cc
      )
endif()

instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
    CRKSPHCorrectionParams.hh
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHSpecificThermalEnergyPolicy.hh
    CRKSPHUtilities.hh
    CRKSPHUtilitiesInline.hh
    ContinuityVolumePolicy.hh
    HullVolumePolicy.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    SolidCRKSPHHydroBaseRZInline.hh
    SurfaceNodeCoupling.hh
    centerOfMass.hh
    computeCRKSPHCorrections.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHMoments.hh
    computeCRKSPHSumMassDensity.hh
    computeCRKSPHSumVolume.hh
    computeHVolumes.hh
    computeHullSumMassDensity.hh
    computeHullVolumes.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    computeVoronoiVolume.hh
    detectSurface.hh
    #flagSurfaceNeighbors.hh
    gradientCRKSPH.hh
    interpolateCRKSPH.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources}
                HEADERS     ${CRKSPH_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_CRKSPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )
      
list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CRKSPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)