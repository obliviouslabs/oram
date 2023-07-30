get_property(INCLUDE_FILES GLOBAL PROPERTY INCLUDE_FILES_G)
get_property(SOURCE_FILES_NM GLOBAL PROPERTY SOURCE_FILES_NM_G)

add_custom_target(
  clangformat
  COMMAND /usr/bin/clang-format
  -style=LLVM
  -i
  ${SOURCE_FILES_NM}
)