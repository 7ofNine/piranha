# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# CMake module to setup CPack for piranha.

# Let's build a tarball under *NIX, and an auto-installing NSIS
# package elsewhere (i.e., Windows).
IF(UNIX)
	SET(CPACK_GENERATOR "TGZ")
ELSE(UNIX)
	SET(CPACK_GENERATOR "NSIS")
ENDIF(UNIX)

SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(CPACK_PACKAGE_VERSION "${PIRANHA_VERSION_STRING}")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Piranha - A Computer Algebra System for Celestial Mechanics")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_STRIP_FILES TRUE)
#SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${PIRANHA_VERSION_STRING}-${CPACK_SYSTEM_NAME}")

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
		CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Examples.lnk\\\" \\\"$INSTDIR\\\\examples\\\"
		CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\License.lnk\\\" \\\"$INSTDIR\\\\license.txt\\\"
		CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Changelog.lnk\\\" \\\"$INSTDIR\\\\changelog.txt\\\"
		"
	)
	# Delete shortcuts when uninstalling.
	SET(CPACK_NSIS_DELETE_ICONS_EXTRA
		"
		Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Pyranha.lnk\\\"
		Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Examples.lnk\\\"
		Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\License.lnk\\\"
		Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Changelog.lnk\\\"
		"
	)
ENDIF(CPACK_GENERATOR MATCHES "NSIS")
