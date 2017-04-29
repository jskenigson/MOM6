#

SHELL = bash
COMPILERS = gnu intel pgi

FMS_SRC = MOM6-examples/src/FMS
MOM6_SRC = MOM6-examples/src/MOM6
SIS2_SRC = MOM6-examples/src/SIS2
ICEBERGS_SRC = MOM6-examples/src/icebergs
MKMF_SRC = MOM6-examples/src/mkmf
LIST_PATHS = $(MKMF_SRC)/bin/list_paths
MKMF = $(MKMF_SRC)/bin/mkmf
BUILD = build
SITE = ncrc

POSSIBLE_DYNAMIC_CONFIGURATIONS=ocean_only ice_ocean_SIS ice_ocean_SIS2 land_ice_ocean_LM3_SIS2 coupled_AM2_LM3_SIS coupled_AM2_LM3_SIS2
#CPP_DEFS = -D_FILE_VERSION="`../../../../../$(MKMF_SRC)/bin/git-version-string $$<`"

# Converts a path a/b/c to a list "a b c"
slash_to_list = $(subst /, ,$(1))
# Replaces a path a/b/c with ../../../
noop =
rel_path = $(subst $(noop) $(noop),,$(patsubst %,../,$(call slash_to_list,$(1))))
# Returns REPRO=1, DEBUG=1 or COVERAGE=1 for repro, debug or coverage, respectively
make_args = $(subst repro,REPRO,$(subst debug,DEBUG,$(subst coverage,COVERAGE,$(1)=1)))

# Environments (env files are source before compilation)
define make-env-dep
$(BUILD)/$(1)/$(2)/fms/libfms.a: $(BUILD)/$(c)/env
endef
$(foreach c, $(COMPILERS),$(foreach m,repro debug coverage,$(eval $(call make-env-dep,$c,$m))))
ifeq ($(SITE),ncrc)
$(BUILD)/gnu/env:
	mkdir -p $(@D)
	echo 'module unload PrgEnv-intel netcdf ; module load gcc/4.9.3 PrgEnv-gnu cray-netcdf' > $@
endif

# build/compiler/mode/fms/libfms.a
# fms_compiler = gnu, intel, pgi, cray, ...
# fms_mode = repro, debug, coverage, ...
fms_compiler = $(word 2,$(call slash_to_list, $(1)))
fms_mode = $(word 3,$(call slash_to_list, $(1)))
$(BUILD)/%/fms/path_names $(BUILD)/%/fms/Makefile $(BUILD)/%/fms/libfms.a: $(wildcard $(FMS_SRC)/*) $(wildcard $(FMS_SRC)/*/*) $(wildcard $(FMS_SRC)/*/*/*)
	mkdir -p $(@D)
	(cd $(@D); rm -f path_names; $(call rel_path,$(@D))$(LIST_PATHS) $(call rel_path,$(@D))$(FMS_SRC))
	(cd $(@D); $(call rel_path,$(@D))$(MKMF) -t $(call rel_path,$(@D))$(MKMF_SRC)/templates/$(SITE)-$(call fms_compiler,$@).mk -p libfms.a -c '-Duse_netCDF' path_names)
	(cd $(@D); source ../../env && make NETCDF=3 $(call make_args, $(call fms_mode, $@)) libfms.a)

# build/compiler/mode/icebergs/libicebergs.a
$(BUILD)/%/icebergs/path_names $(BUILD)/%/icebergs/Makefile $(BUILD)/%/icebergs/libicebergs.a: $(wildcard $(ICEBERGS_SRC)/*) $(wildcard $(ICEBERGS_SRC)/*/*) $(wildcard $(ICEBERGS_SRC)/*/*/*)
	mkdir -p $(@D)
	(cd $(@D); rm -f path_names; $(call rel_path,$(@D))$(LIST_PATHS) $(call rel_path,$(@D))$(ICEBERGS_SRC))
	(cd $(@D); $(call rel_path,$(@D))$(MKMF) -t $(call rel_path,$(@D))$(MKMF_SRC)/templates/$(SITE)-$(call fms_compiler,$@).mk -p libicebergs.a -o '-I../fms/$(EXEC_MODE)' path_names)
	(cd $(@D); source ../../env && make NETCDF=3 $(call make_args, $(call fms_mode, $@)) libicebergs.a)


# Generate list of MOM6 executables that depend on libfms.a
define make-fms-dep
$(BUILD)/$(1)/$(2)/libmom6.a $(BUILD)/$(1)/$(2)/MOM6: $(BUILD)/$(1)/fms/libfms.a
endef
$(foreach c,$(COMPILERS),$(foreach m,repro debug coverage,$(foreach d,dynamic dynamic_symmetric,$(foreach o,$(POSSIBLE_DYNAMIC_CONFIGURATIONS) libmom,$(eval $(call make-fms-dep,$(c)/$(m),$(d)/$(o)))))))
# Generate list of MOM6 executables with specific settings for LIST_PATH_ARGS
define make-memory-dep
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libmom6.a $(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATH_ARGS += $(MOM6_SRC)/config_src/$(3)/ $(MOM6_SRC)/config_src/$(if $(findstring ocean_only,$(4)),solo_driver,coupled_driver)/
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libmom6.a $(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: CPP_DEFS += -DSTATSLABEL=\"$(STATS_PLATFORM)$(1)$(STATS_COMPILER_VER)\"
ifeq ($(4),ice_ocean_SIS2)
$(BUILD)/$(1)/$(2)/$(3)/$(4)/libsis2.a $(BUILD)/$(1)/$(2)/$(3)/$(4)/MOM6: LIST_PATH_ARGS += $(SIS2_SRC)/
endif
endef
$(foreach c,$(COMPILERS),$(foreach m,repro debug coverage,$(foreach d,dynamic dynamic_symmetric,$(foreach o,$(POSSIBLE_DYNAMIC_CONFIGURATIONS) libmom,$(eval $(call make-memory-dep,$(c),$(m),$(d),$(o)))))))

# build/compiler/mode/memory/configuration/MOM6
# compiler = gnu, intel, pgi, cray, ...
# mode = repro, debug, coverage, ...
# memory = dynamic, dynamic_symmetric, static
# configuration = ocean_only, ice_ocean_SIS, ice_ocean_SIS2, land_ice_ocean_LM3_SIS2, coupled_AM2_LM3_SIS, coupled_AM2_LM3_SIS2
# configuration = CVmix_SCM_tests DOME Phillips_2layer SCM_idealized_hurricane adjustment2d benchmark buoy_forced_basin circle_obcs double_gyre ...
compiler = $(word 2,$(call slash_to_list, $(1)))
mode = $(word 3,$(call slash_to_list, $(1)))
memory = $(word 4,$(call slash_to_list, $(1)))
configuration = $(word 5,$(call slash_to_list, $(1)))
LIST_PATH_ARGS = $(MOM6_SRC)/src/*/ $(MOM6_SRC)/src/*/*/
$(BUILD)/%/path_names $(BUILD)/%/Makefile $(BUILD)/%/MOM6 $(BUILD)/%/libmom6.a: $(MOM6_CONFIG_SRC) $(wildcard $(foreach d,$(LIST_PATH_ARGS),$(d)*.F90 $(d)*.h))
	mkdir -p $(@D)
	(cd $(@D); rm -f path_names; $(call rel_path,$(@D))$(LIST_PATHS) $(foreach p,$(LIST_PATH_ARGS),$(call rel_path,$(@D))$(p)))
	(cd $(@D); $(call rel_path,$(@D))$(MKMF) -t $(call rel_path,$(@D))$(MKMF_SRC)/templates/$(SITE)-$(call compiler,$@).mk -p $(@F) -o '-I../../fms/$(EXEC_MODE)' -l '-L../../fms/$(EXEC_MODE) -lfms' -c '$(CPP_DEFS)' path_names)
	(cd $(@D); source ../../../env && make $(call make_args, $(call fms_mode, $@)) $(@F))

MOM6-examples:
	git clone --recursive https://github.com/NOAA-GFDL/MOM6-examples.git
