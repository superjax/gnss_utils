
get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/gnss_utils-targets.cmake)
get_filename_component(aysnc_comm_INCLUDE_DIRS "${SELF_DIR}/../../include/gnss_utils" ABSOLUTE)
set(gnss_utils_LIBRARIES gnss_utils)