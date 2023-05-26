if [ "$EUID" -eq 0 ]
  then echo "This script is not intended to be run as user root"
  exit 1
fi

EXTENSION_TARGET="${HOME}/klipper/klippy/extras"
declare -a SCRIPTS=("linear_movement_vibrations.py" "linear_movement_plot_lib_stat.py")
GIT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo -e "Installation Script for Klipper Linear Movement Analysis by MarschallMarc#6420"
echo -e "Linking extension file"

for i in "${SCRIPTS[@]}"
do
    if [ -d "${EXTENSION_TARGET}" ]; then
        rm -f "${EXTENSION_TARGET}/$i"
        ln -sf "${GIT_DIR}/$i" "${EXTENSION_TARGET}/$i"
    else
        echo -e "${EXTENSION_TARGET} not found, exiting installation."
        exit 1
    fi
done

~/klippy-env/bin/pip install -v matplotlib scipy