MFEM_DIR ?= $(APPLICATION_DIR)/../mfem/installed

include $(MFEM_DIR)/share/mfem/config.mk

NEW_MFEM_INCFLAGS = $(shell echo "$(MFEM_INCFLAGS)" | sed -E 's/COMPILE_LANGUAGE:CXX>:SHELL:[^>]*>//g')

ADDITIONAL_INCLUDES += $(NEW_MFEM_INCFLAGS)
ADDITIONAL_LIBS     += $(MFEM_LIBS) -lmfem-common