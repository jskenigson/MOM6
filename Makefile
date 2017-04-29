#

SHELL = bash
COMPILERS = gnu intel pgi

SRC_DIR = MOM6-examples/src
FMS_SRC = $(SRC_DIR)/FMS
MOM6_SRC = $(SRC_DIR)/MOM6
SIS2_SRC = $(SRC_DIR)/SIS2
ICE_OCEAN_EXTRAS_SRC = $(SRC_DIR)/ice_ocean_extras
ICEBERGS_SRC = $(SRC_DIR)/icebergs
MKMF_SRC = $(SRC_DIR)/mkmf
LIST_PATHS = $(MKMF_SRC)/bin/list_paths
MKMF = $(MKMF_SRC)/bin/mkmf
BUILD = build
SITE = ncrc

POSSIBLE_DYNAMIC_CONFIGURATIONS=ocean_only ice_ocean_SIS ice_ocean_SIS2 land_ice_ocean_LM3_SIS2 coupled_AM2_LM3_SIS coupled_AM2_LM3_SIS2
POSSIBLE_SIS2_CONFIGURATIONS=ice_ocean_SIS2 coupled_AM2_LM3_SIS2
CPP_DEFS = -D_FILE_VERSION="`../../../../../$(MKMF_SRC)/bin/git-version-string $$<`"

# Converts a path a/b/c to a list "a b c"
slash_to_list = $(subst /, ,$(1))
# Replaces a path a/b/c with ../../../
noop =
rel_path = $(subst $(noop) $(noop),,$(patsubst %,../,$(call slash_to_list,$(1))))
# Returns REPRO=1, DEBUG=1 or COVERAGE=1 for repro, debug or coverage, respectively
make_args = $(subst repro,REPRO,$(subst debug,DEBUG,$(subst coverage,COVERAGE,$(1)=1)))

# Environments (env files are source before compilation)
# make-env-dep: $(1) = compiler, $(2) = mode
define make-env-dep
$(BUILD)/$(1)/$(2)/fms/libfms.a: $(BUILD)/$(1)/env
endef
$(foreach c, $(COMPILERS),$(foreach m,repro debug coverage,$(eval $(call make-env-dep,$c,$m))))
ifeq ($(SITE),ncrc)
$(BUILD)/gnu/env:
	mkdir -p $(@D)
	echo 'module unload PrgEnv-intel netcdf ; module load gcc/4.9.3 PrgEnv-gnu cray-netcdf' > $@
endif

# path_names:
# - must have LIST_PATH_ARGS set for final target
$(BUILD)/%/path_names: $(wildcard $(LIST_PATHS_ARGS)*) $(wildcard $(LIST_PATHS_ARGS)*/*) $(wildcard $(LIST_PATHS_ARGS)*/*/*) $(wildcard $(LIST_PATHS_ARGS)*/*/*/*)
$(BUILD)/%/path_names:
	mkdir -p $(@D)
	(cd $(@D); rm -f path_names; $(call rel_path,$(@D))$(LIST_PATHS) $(foreach p,$(LIST_PATHS_ARGS),$(call rel_path,$(@D))$(p)))

# Makefile:
# - must have MKMF_OPTS set for final target
# fms_compiler = gnu, intel, pgi, cray, ...
fms_compiler = $(word 2,$(call slash_to_list, $(1)))
$(BUILD)/%/Makefile: $(BUILD)/%/path_names
	(cd $(@D); $(call rel_path,$(@D))$(MKMF) -t $(call rel_path,$(@D))$(MKMF_SRC)/templates/$(SITE)-$(call fms_compiler,$@).mk $(MKMF_OPTS) -c '-Duse_netCDF' path_names)

# build/compiler/mode/fms/libfms.a
# fms_mode = repro, debug, coverage, ...
fms_mode = $(word 3,$(call slash_to_list, $(1)))
$(BUILD)/%/fms/path_names: LIST_PATHS_ARGS = $(FMS_SRC)/
$(BUILD)/%/fms/Makefile: MKMF_OPTS = -p libfms.a
$(BUILD)/%/libfms.a: $(BUILD)/%/Makefile
	(cd $(@D); source ../../env && make NETCDF=3 $(call make_args, $(call fms_mode, $@)) $(@F))

# build/compiler/mode/icebergs/libicebergs.a
$(BUILD)/%/icebergs/path_names: LIST_PATHS_ARGS = $(ICEBERGS_SRC)/
$(BUILD)/%/icebergs/Makefile: MKMF_OPTS = -p libicebergs.a -o '-I../fms'
$(BUILD)/%/icebergs/libicebergs.a: $(BUILD)/%/icebergs/Makefile $(BUILD)/%/fms/libfms.a
	(cd $(@D); source ../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))

# build/compiler/mode/ice_ocean_extras/libice_ocean_extras.a
$(BUILD)/%/ice_ocean_extras/path_names: LIST_PATHS_ARGS = $(ICE_OCEAN_EXTRAS_SRC)/ $(FMS_SRC)/include/fms_platform.h
$(BUILD)/%/ice_ocean_extras/Makefile: MKMF_OPTS = -p libice_ocean_extras.a -o '-I../fms'
$(BUILD)/%/ice_ocean_extras/libice_ocean_extras.a: $(BUILD)/%/ice_ocean_extras/Makefile $(BUILD)/%/fms/libfms.a
	(cd $(@D); source ../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))

# build/compiler/mode/sis2/libsis2.a
$(BUILD)/%/sis2/path_names: LIST_PATHS_ARGS += $(SIS2_SRC)/
$(BUILD)/%/sis2/Makefile: MKMF_OPTS = -p libsis2.a -o '-I../../fms -I../mom6 -I../../icebergs -I../../ice_ocean_extras'
$(BUILD)/%/sis2/libsis2.a: $(BUILD)/%/sis2/Makefile $(BUILD)/%/mom6/libmom6.a
	(cd $(@D); source ../../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))

# Generate list of MOM6 executables that depend on libfms.a
# make-fms-dep: $(1) = compiler + mode, $(2) = memory style + mom6 configuration
define make-fms-dep
$(BUILD)/$(1)/$(2)/libmom6.a $(BUILD)/$(1)/$(2)/MOM6: $(BUILD)/$(1)/fms/libfms.a
endef
$(foreach c,$(COMPILERS),$(foreach m,repro debug coverage,$(foreach d,dynamic dynamic_symmetric,$(foreach o,$(POSSIBLE_DYNAMIC_CONFIGURATIONS) mom6,$(eval $(call make-fms-dep,$(c)/$(m),$(d)/$(o)))))))
# Generate lists of LIST_PATH_ARGS for MOM6 executables and libraries
# $(1) = compiler, $(2) = mode, $(3) = memory style, $(4) = mom6 configuration
#$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATHS_ARGS += $(MOM6_SRC)/config_src/$(3)/ $(MOM6_SRC)/config_src/$(if $(findstring ocean_only,$(4)),solo_driver,coupled_driver)/
define mom6-variables
ifeq ($(4),mom6)
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libmom6.a: LIST_PATHS_ARGS += $(MOM6_SRC)/src/*/ $(MOM6_SRC)/src/*/*/
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libmom6.a: LIST_PATHS_ARGS += $(MOM6_SRC)/config_src/$(3)/ $(MOM6_SRC)/config_src/coupled_driver/
endif
ifeq ($(4),ocean_only)
$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATHS_ARGS += $(MOM6_SRC)/src/*/ $(MOM6_SRC)/src/*/*/
$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATHS_ARGS += $(MOM6_SRC)/config_src/$(3)/ $(MOM6_SRC)/config_src/solo_driver/
endif
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libmom6.a $(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: CPP_DEFS += -DSTATSLABEL=\"$(STATS_PLATFORM)$(1)$(STATS_COMPILER_VER)\"
endef
$(foreach c,$(COMPILERS),$(foreach m,repro debug coverage,$(foreach d,dynamic dynamic_symmetric,$(foreach o,mom6 $(POSSIBLE_DYNAMIC_CONFIGURATIONS),$(eval $(call mom6-variables,$(c),$(m),$(d),$(o)))))))
# Generate lists of LIST_PATH_ARGS for SIS2 executables and libraries
# $(1) = compiler, $(2) = mode, $(3) = memory style, $(4) = mom6 configuration
define sis2-variables
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libsis2.a: $(BUILD)/$(1)/$(2)/icebergs/libicebergs.a $(BUILD)/$(1)/$(2)/ice_ocean_extras/libice_ocean_extras.a
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libsis2.a: LIST_PATHS_ARGS += $(MOM6_SRC)/src/framework/MOM_memory_macros.h
ifeq ($(4),ice_ocean_SIS2)
$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATHS_ARGS += $(FMS_SRC)/coupler/
$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: MKMF_OPTS += XXXXX
endif
endef
#$(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATHS_ARGS += $(SIS2_SRC)/
$(foreach c,$(COMPILERS),$(foreach m,repro debug coverage,$(foreach d,dynamic dynamic_symmetric,$(foreach o,sis2 $(POSSIBLE_SIS2_CONFIGURATIONS),$(eval $(call sis2-variables,$(c),$(m),$(d),$(o)))))))

# build/compiler/mode/mom6_memory/mom6_configuration/MOM6
# compiler = gnu, intel, pgi, cray, ...
# mode = repro, debug, coverage, ...
# mom6_memory = dynamic, dynamic_symmetric, static
# mom6_configuration = ocean_only, ice_ocean_SIS, ice_ocean_SIS2, land_ice_ocean_LM3_SIS2, coupled_AM2_LM3_SIS, coupled_AM2_LM3_SIS2
# mom6_configuration = CVmix_SCM_tests DOME Phillips_2layer SCM_idealized_hurricane adjustment2d benchmark buoy_forced_basin circle_obcs double_gyre ...
#compiler = $(word 2,$(call slash_to_list, $(1)))
#mode = $(word 3,$(call slash_to_list, $(1)))
#mom6_memory = $(word 4,$(call slash_to_list, $(1)))
#mom6_configuration = $(word 5,$(call slash_to_list, $(1)))
$(BUILD)/%/mom6/libmom6.a: MKMF_OPTS = -p libmom6.a -o '-I../../fms' -l '-L../../fms -lfms' -c '$(CPP_DEFS)'
$(BUILD)/%/mom6/libmom6.a: $(BUILD)/%/mom6/Makefile
	(cd $(@D); source ../../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))
$(BUILD)/%/MOM6: MKMF_OPTS += -p MOM6 -o '-I../../fms' -l '-L../../fms -lfms' -c '$(CPP_DEFS)'
$(BUILD)/%/MOM6: $(BUILD)/%/Makefile
	(cd $(@D); source ../../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))

MOM6-examples:
	git clone --recursive https://github.com/NOAA-GFDL/MOM6-examples.git

whats_built:
	find $(BUILD) -name "MOM6" -o -name "lib*.a"

test: \
build/gnu/repro/dynamic/ocean_only/MOM6 \
build/gnu/repro/dynamic/mom6/libmom6.a \
build/gnu/repro/dynamic_symmetric/ocean_only/MOM6 \
build/gnu/repro/dynamic_symmetric/mom6/libmom6.a
