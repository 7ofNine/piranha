#! /bin/sh

CURRENT_RELEASE=`git log --pretty=oneline HEAD^1..HEAD | head -n 1 | awk '{print $1}'`
PREVIOUS_RELEASE=`cat current_release`
cp Changelog Changelog.tmp
echo "`date +%Y.%m.%d` ${CURRENT_RELEASE}" > Changelog
echo "" >> Changelog
git log --pretty=oneline ${PREVIOUS_RELEASE}.. | grep -i "clog:" | awk -F "clog:" '{print $2}' | sed 's/^[ \t]*//;s/[ \t]*$//'|sed 's/^/\*\ /' >> Changelog
echo "" >> Changelog
cat Changelog.tmp >> Changelog
echo ${CURRENT_RELEASE} > current_release
rm Changelog.tmp
