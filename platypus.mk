MFEM_DIR				:=$(APPLICATION_DIR)/../mfem/build

include $(MFEM_DIR)/config/config.mk
MFEM_INCLUDES 		:= -I$(MFEM_INC_DIR)/config -I$(MFEM_DIR)/ -I$(MFEM_DIR)/../miniapps/common
MFEM_LIBS 			:= -L$(MFEM_DIR) -lmfem -lrt -L$(MFEM_DIR)/miniapps/common -lmfem-common $(MFEM_LIB)

PLATYPUS_LIB		:= $(PLATYPUS_DIR)/lib

SPDLOG_LIB 			 :=$(PLATYPUS_DIR)/build/_deps/spdlog-build/
SPDLOG_INCLUDE 		 :=$(PLATYPUS_DIR)/build/_deps/spdlog-src/include/

PLATYPUS_INCLUDE_LIB := $(sort $(dir $(shell find $(PLATYPUS_DIR)/src/ -name "*.hpp")))
PLATYPUS_INCLUDE	 := $(foreach d, $(PLATYPUS_INCLUDE_LIB), -I$d) -I$(SPDLOG_INCLUDE)

PLATYPUS_CXX_FLAGS 	:= -Wall $(MFEM_INCLUDES) -I$(MFEM_INC_DIR)/config
PLATYPUS_LDFLAGS	:= -Wl,-rpath, $(MFEM_LIBS) -L$(SPDLOG_LIB) -lspdlog

libmesh_CXXFlags 		+= $(PLATYPUS_CXX_FLAGS)
libmesh_LDFLAGS 		+= $(PLATYPUS_LDFLAGS)

ADDITIONAL_INCLUDES 	+= $(PLATYPUS_CXX_FLAGS)
ADDITIONAL_LIBS 		+= $(PLATYPUS_LDFLAGS)

$(info ADDITIONAL_INCLUDES = $(ADDITIONAL_INCLUDES));
$(info ADDITIONAL_LIBS     = $(ADDITIONAL_LIBS));