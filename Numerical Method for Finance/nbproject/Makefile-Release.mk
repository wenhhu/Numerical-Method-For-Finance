#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=clang
CCC=clang++
CXX=clang++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=CLang-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/57789b2f/Binomialpricing.o \
	${OBJECTDIR}/_ext/57789b2f/BlackScholes.o \
	${OBJECTDIR}/_ext/57789b2f/Cholesky.o \
	${OBJECTDIR}/_ext/57789b2f/Functions.o \
	${OBJECTDIR}/_ext/57789b2f/NumMethod.o \
	${OBJECTDIR}/_ext/57789b2f/iterative_linear_solver.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/numerical_method_for_finance

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/numerical_method_for_finance: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/numerical_method_for_finance ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/57789b2f/Binomialpricing.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Binomialpricing.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/Binomialpricing.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Binomialpricing.cpp

${OBJECTDIR}/_ext/57789b2f/BlackScholes.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/BlackScholes.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/BlackScholes.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/BlackScholes.cpp

${OBJECTDIR}/_ext/57789b2f/Cholesky.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Cholesky.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/Cholesky.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Cholesky.cpp

${OBJECTDIR}/_ext/57789b2f/Functions.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Functions.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/Functions.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/Functions.cpp

${OBJECTDIR}/_ext/57789b2f/NumMethod.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/NumMethod.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/NumMethod.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/NumMethod.cpp

${OBJECTDIR}/_ext/57789b2f/iterative_linear_solver.o: /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/iterative_linear_solver.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/57789b2f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/57789b2f/iterative_linear_solver.o /Users/wenhaohu/NetBeansProjects/Numerical_Method_For_Finace/Numerical\ Method\ for\ Finance/iterative_linear_solver.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/numerical_method_for_finance

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
