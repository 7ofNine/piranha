#! /bin/bash

if [[ "$@" == "" ]]
then for i in `find ./src -regex "\(.*\.h\\|.*\.cpp\)"`
    do bcpp -i 2 "${i}" > "${i}".new
    mv "${i}".new "${i}"
  done
else
  for i in "$@"
    do bcpp -i 2 "${i}" > "${i}".new
    mv "${i}".new "${i}"
  done
fi
