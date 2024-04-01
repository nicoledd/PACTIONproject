# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_ALL_INSTALL_TYPES "Full;Developer")
set(CPACK_BUILD_SOURCE_DIRS "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1;/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENTS_ALL "headers;library;html_documentation;bin")
set(CPACK_COMPONENTS_ALL_SET_BY_USER "TRUE")
set(CPACK_COMPONENT_BIN_DESCRIPTION "Command line utilities")
set(CPACK_COMPONENT_BIN_DISPLAY_NAME "Command line utilities")
set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "Components needed to develop software using LEMON")
set(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION "Documentation of LEMON")
set(CPACK_COMPONENT_HEADERS_DEPENDS "library")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ header files")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ headers")
set(CPACK_COMPONENT_HEADERS_GROUP "Development")
set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES "Developer;Full")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION "Doxygen generated documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES "Full")
set(CPACK_COMPONENT_LIBRARY_DESCRIPTION "DLL and import library")
set(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME "Dynamic-link library")
set(CPACK_COMPONENT_LIBRARY_GROUP "Development")
set(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES "Developer;Full")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/opt/homebrew/Cellar/cmake/3.27.4/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "LEMON built using CMake")
set(CPACK_DMG_SLA_USE_RESOURCE_FILE_LICENSE "ON")
set(CPACK_GENERATOR "NSIS")
set(CPACK_INNOSETUP_ARCHITECTURE "x64")
set(CPACK_INSTALL_CMAKE_PROJECTS "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build;LEMON;ALL;/")
set(CPACK_INSTALL_PREFIX "/Users/nsdong2/lemon")
set(CPACK_MODULE_PATH "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/cmake")
set(CPACK_NSIS_CONTACT "lemon-user@lemon.cs.elte.hu")
set(CPACK_NSIS_CREATE_ICONS_EXTRA "
    CreateShortCut \"$SMPROGRAMS\\$STARTMENU_FOLDER\\Documentation.lnk\" \"$INSTDIR\\share\\doc\\index.html\"
    ")
set(CPACK_NSIS_DELETE_ICONS_EXTRA "
    !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    Delete \"$SMPROGRAMS\\$MUI_TEMP\\Documentation.lnk\"
    ")
set(CPACK_NSIS_DISPLAY_NAME "LEMON 1.3.1 LEMON")
set(CPACK_NSIS_DISPLAY_NAME_SET "TRUE")
set(CPACK_NSIS_HELP_LINK "http:\\\\lemon.cs.elte.hu")
set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\lemon.ico")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_MUI_ICON "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/cmake/nsis/lemon.ico")
set(CPACK_NSIS_MUI_UNIICON "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/cmake/nsis/uninstall.ico")
set(CPACK_NSIS_PACKAGE_NAME "LEMON 1.3.1 LEMON")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\lemon.cs.elte.hu")
set(CPACK_OBJDUMP_EXECUTABLE "/Library/Developer/CommandLineTools/usr/bin/objdump")
set(CPACK_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk")
set(CPACK_OUTPUT_CONFIG_FILE "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/opt/homebrew/Cellar/cmake/3.27.4/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LEMON - Library for Efficient Modeling and Optimization in Networks")
set(CPACK_PACKAGE_FILE_NAME "LEMON-1.3.1-Darwin")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "LEMON 1.3.1")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "LEMON 1.3.1")
set(CPACK_PACKAGE_NAME "LEMON")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "EGRES")
set(CPACK_PACKAGE_VERSION "1.3.1")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/LICENSE")
set(CPACK_RESOURCE_FILE_README "/opt/homebrew/Cellar/cmake/3.27.4/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/opt/homebrew/Cellar/cmake/3.27.4/share/cmake/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TBZ2;TGZ;TXZ;TZ")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/CPackSourceConfig.cmake")
set(CPACK_SOURCE_RPM "OFF")
set(CPACK_SOURCE_TBZ2 "ON")
set(CPACK_SOURCE_TGZ "ON")
set(CPACK_SOURCE_TXZ "ON")
set(CPACK_SOURCE_TZ "ON")
set(CPACK_SOURCE_ZIP "OFF")
set(CPACK_SYSTEM_NAME "Darwin")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Darwin")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
