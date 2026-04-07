#ifndef VTU_INCLUDE
#define VTU_INCLUDE(a,b,c)
#endif


VTU_INCLUDE("position",Y%pos,3)
VTU_INCLUDE("velocity",Y%vel,3)
VTU_INCLUDE("Temp",Y%T,1)
VTU_INCLUDE("AngularVelocity",Y%ang_vel,3)
VTU_INCLUDE("Metal",Y%metal,1)
VTU_INCLUDE("Oxide",Y%oxide,1)

VTU_INCLUDE("RHOp",rprop%RHOp,1)
VTU_INCLUDE("Dp",rprop%Dp,1)
VTU_INCLUDE("VOLp",rprop%VOLp,1)
VTU_INCLUDE("JDPi",rprop%JDPi,1)
VTU_INCLUDE("JDPe",rprop%JDPe,1)
VTU_INCLUDE("JSPL",rprop%JSPL,1)
VTU_INCLUDE("JSPT",rprop%JSPT,1)
VTU_INCLUDE("FLUCTF",rprop%FLUCTF,3)
VTU_INCLUDE("WDOT",rprop%WDOT,3)
VTU_INCLUDE("IDp",rprop%IDp,1)
VTU_INCLUDE("BRNT",rprop%BRNT,1)
VTU_INCLUDE("XIPAR",rprop%XIPAR,1)
VTU_INCLUDE("XIPERP",rprop%XIPERP,1)
VTU_INCLUDE("XIT",rprop%XIT,1)


VTU_INCLUDE("tag",tag,3)
VTU_INCLUDE("ParticleRank",iprop%ParticleRank,1)
VTU_INCLUDE("xBin",iprop%xBin,1)
VTU_INCLUDE("yBin",iprop%yBin,1)
VTU_INCLUDE("zBin",iprop%zBin,1)
VTU_INCLUDE("RemoveParticle",iprop%RemoveParticle,1)
