# CMake module to setup CPack for piranha.

# Let's build a tarball under *NIX, and an auto-installing NSIS
# package elsewhere (i.e., Windows).
IF(UNIX)
	SET(CPACK_GENERATOR "TGZ")
ELSE(UNIX)
	SET(CPACK_GENERATOR "NSIS")
ENDIF(UNIX)

SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")

IF(CPACK_GENERATOR MATCHES "NSIS")
	# With NSIS we have further possibilities for customization.
	#SET(CPACK_NSIS_DISPLAY_NAME "${PROJECT_NAME} ${VERSION}")
	# Apparently this escaping madness is necessary due to an NSIS bug.
	SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\piranha.tuxfamily.org")
	SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\piranha.tuxfamily.org")
	# Add shortcuts to the Start Menu.
	SET(CPACK_NSIS_CREATE_ICONS_EXTRA
		"
		CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Pyranha.lnk\\\" \\\"$INSTDIR\\\\Console.exe\\\"
		"
	)
	# Delete shortcuts when uninstalling.
	SET(CPACK_NSIS_DELETE_ICONS_EXTRA
		"
		Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Pyranha.lnk\\\"
		"
	)
ENDIF(CPACK_GENERATOR MATCHES "NSIS")
