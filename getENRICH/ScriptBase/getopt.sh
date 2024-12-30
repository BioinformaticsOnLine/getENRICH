#!/bin/bash

#Store the variable name
function getVName {
    local OPT_DESC=${1}
    local OPTION=${2}
    local VAR=$(echo ${OPT_DESC} | sed -e "s/.*\[\?-${OPTION} \([A-Z_]\+\).*/\1/g" -e "s/.*\[\?-\(${OPTION}\).*/\1FLAG/g")

    if [[ "${VAR}" == "${1}" ]]; then
        echo ""
    else
        echo ${VAR}
    fi
}
#Parsing the options
function parse_options {
    local OPT_DESC=${1}
    local INPUT=$(getInOpt "${OPT_DESC}")

    shift
    while getopts ${INPUT} OPTION ${@};
    do
        [ ${OPTION} == "?" ] && usage
        VARNAME=$(getVName "${OPT_DESC}" "${OPTION}")
            [ "${VARNAME}" != "" ] && eval "${VARNAME}=${OPTARG:-true}" 
    done

    checkR "${OPT_DESC}"

}
#Checking the options
function checkR {
    local OPT_DESC=${1}
    local REQUIRED=$(getR "${OPT_DESC}" | sed -e "s/\://g")
    while test -n "${REQUIRED}"; do
        OPTION=${REQUIRED:0:1}
        VARNAME=$(getVName "${OPT_DESC}" "${OPTION}")
                [ -z "${!VARNAME}" ] && printf "ERROR: %s\n" "Option -${OPTION} must been set." && usage
        REQUIRED=${REQUIRED:1}
    done
}

#Input option dealing 
function getInOpt {
    local OPT_DESC=${1}
    echo ${OPT_DESC} | sed -e "s/\([a-zA-Z]\) [A-Z_]\+/\1:/g" -e "s/[][ -]//g"
}

#Option dealing
function getO {
    local OPT_DESC=${1}
    echo ${OPT_DESC} | sed -e "s/[^[]*\(\[[^]]*\]\)[^[]*/\1/g" -e "s/\([a-zA-Z]\) [A-Z_]\+/\1:/g" -e "s/[][ -]//g"
}

#Required dealing
function getR {
    local OPT_DESC=${1}
    echo ${OPT_DESC} | sed -e "s/\([a-zA-Z]\) [A-Z_]\+/\1:/g" -e "s/\[[^[]*\]//g" -e "s/[][ -]//g"
}

#print the usage
function usage {
    printf "Usage:\n\t%s\n" "${0} ${OPT_DESC}"
    exit 10
}

