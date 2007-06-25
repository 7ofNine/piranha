#! /bin/bash

for i in `find ./ -regex "\(.*\.h\\|.*\.cpp\)"`
  do bcpp -i 2 "${i}" > "${i}".new
  mv "${i}".new "${i}"
done
