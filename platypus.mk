MFEM_DIR			:=$(APPLICATION_DIR)/../mfem/build

include $(MFEM_DIR)/config/config.mk
MFEM_INCLUDES 		:= -I$(MFEM_INC_DIR)/config -I$(MFEM_DIR)/ -I$(MFEM_DIR)/../miniapps/common $(MFEM_INCFLAGS)
MFEM_LIBS 			:= -L$(MFEM_DIR) -lmfem -lrt -L$(MFEM_DIR)/miniapps/common -lmfem-common $(MFEM_LIB)

ADDITIONAL_INCLUDES += $(MFEM_INCLUDES)
ADDITIONAL_LIBS 	+= -Wl, $(MFEM_LIBS) ${MFEM_EXT_LIBS}

$(info ADDITIONAL_INCLUDES = $(ADDITIONAL_INCLUDES));
$(info ADDITIONAL_LIBS     = $(ADDITIONAL_LIBS));