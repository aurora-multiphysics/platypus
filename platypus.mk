MFEM_DIR			?=$(APPLICATION_DIR)/../mfem/build
MFEM_INC_DIR                    ?=$(MFEM_DIR)/include

MFEM_INCLUDES 		:= -I$(MFEM_DIR)/include
MFEM_LIBS 			:= -Wl,-rpath,$(MFEM_DIR)/lib -L$(MFEM_DIR)/lib -lmfem -lrt -lmfem-common

ADDITIONAL_INCLUDES += $(MFEM_INCLUDES)
ADDITIONAL_LIBS 	+= $(MFEM_LIBS) ${MFEM_EXT_LIBS}

$(info ADDITIONAL_INCLUDES = $(ADDITIONAL_INCLUDES));
$(info ADDITIONAL_LIBS     = $(ADDITIONAL_LIBS));
