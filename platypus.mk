MFEM_DIR			:=$(APPLICATION_DIR)/../mfem/build

include $(MFEM_DIR)/config/config.mk

# Remove cmake CXX flags
NEW_MFEM_INCFLAGS = $(shell echo "$(MFEM_INCFLAGS)" | sed -E 's/COMPILE_LANGUAGE:CXX>:SHELL:[^>]*>//g')

MFEM_INCLUDES 		:= -I$(MFEM_INC_DIR)/config -I$(MFEM_DIR)/ -I$(MFEM_DIR)/../miniapps/common $(NEW_MFEM_INCFLAGS)
MFEM_LIBS 			:= -L$(MFEM_DIR) -lmfem -lrt -L$(MFEM_DIR)/miniapps/common -lmfem-common $(MFEM_LIB)

ADDITIONAL_INCLUDES += -Wno-unused-parameter -Wno-unused-variable $(MFEM_INCLUDES)
ADDITIONAL_LIBS 	+= -Wl, $(MFEM_LIBS) ${MFEM_EXT_LIBS}

$(info ADDITIONAL_INCLUDES = $(ADDITIONAL_INCLUDES));
$(info ADDITIONAL_LIBS     = $(ADDITIONAL_LIBS));