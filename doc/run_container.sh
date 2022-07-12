#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -c config_path -m model_path -x x_directory -y y_directory"
   echo -e "\t-c Configuration File Path, look at Config.csv for an example of config formatting, specifies a specific config file. DEFAULT: Config.csv"
   echo -e "\t-m BNGL Model File Path, look at model.bngl for an example, specifies a specific bngl file to read. DEFAULT: model.bngl"
   echo -e "\t-x X data file directory, directory of where data files to be evolved by model are located, make sure to specify a directory not a FILE. DEFAULT: data/X"
   echo -e "\t-y Y data file directory, directory of where observed data files are located, make sure to specify a directory not a FILE. DEFAULT: data/Y"

}

# default values
parameterC="Config.csv"
parameterM="model.bngl"
parameterX="data/X"
parameterY="data/Y"
parameterT="time_steps.csv"
parameterR="true_rates.csv"
# get new values
while getopts "c:m:x:y:" opt
do
   case "$opt" in
      c ) parameterC="$OPTARG" ;;
      m ) parameterM="$OPTARG" ;;
      x ) parameterX="$OPTARG" ;;
      y ) parameterY="$OPTARG" ;;
      t ) parameterT="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      h ) helpFunction ;; # print help if specified
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case some parameters are empty
if [ -z "$parameterC" ] || [ -z "$parameterM" ] || [ -z "$parameterX" ] || [ -z "$parameterY" ] || [ -z "$parameterT" ] || [ -z "$parameterR" ]
then
   echo "Nothing was specified, loading default values!";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Parameters Loaded Below!"
echo "$parameterC"
echo "$parameterM"
echo "$parameterX"
echo "$parameterY"
echo "$parameterT"
echo "$parameterR"

./CyGMM -m $parameterM -c $parameterC -x $parameterX -y $parameterY -t $parameterT -r $parameterR
