#!/bin/bash
if [ "$EUID" -ne 0 ]
  then echo "This script is not intended to be run as user root"
  exit 1
fi

EXTENSION_TARGET="${HOME}/klipper/klippy/extras"
SCRIPT_NAME="linear_movement_vibrations.py"

function link_extension_file{
    if [ -d "${EXTENSION_TARGET}" ]; then
        rm -f "${EXTENSION_TARGET}/${SCRIPT_NAME}"
        ln -sf "${SRCDIR}/klippy_extra/pam.py" "${KLIPPY_EXTRAS}/pam.py"
    else
        echo -e "${EXTENSION_TARGET} not found, exiting installation."
        exit 1
    fi

}


echo -e "Installation Script for Klipper Linear Movement Analysis by MarschallMarc#6420"
echo -e "Linking extension file"
