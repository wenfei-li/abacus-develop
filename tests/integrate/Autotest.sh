#!/bin/bash

# ABACUS executable path
abacus=abacus
# number of MPI processes
np=4
# threshold with unit: eV
threshold=0.0000001
force_threshold=0.0001
stress_threshold=0.001
# check accuracy
ca=8
# specify the test cases file
cases_file=CASES_CPU.txt
# regex of case name
case='^[^#].*_.*$'
# enable AddressSanitizer
sanitize=false

threshold_file="threshold"   
# can specify the threshold for each test case
# threshold file example:
# threshold 0.0000001
# force_threshold 0.0001
# stress_threshold 0.001
# fatal_threshold 1


while getopts a:n:t:c:s:r:f:g flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
        t) threshold=${OPTARG};;
        c) ca=${OPTARG};;
        s) sanitize=${OPTARG};;
        r) case=${OPTARG};;
        f) cases_file=${OPTARG};;
        g) g=true;; #generate test reference
    esac
done

# number of OpenMP threads
nt=$OMP_NUM_THREADS
if [[ -z "$nt" ]]; then
    nt=$(expr `nproc` / ${np})
    export OMP_NUM_THREADS=${nt}
fi

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus"
echo "Number of processes: $np"
echo "Number of threads: $nt"
echo "Test accuracy totenergy: $threshold eV"
echo "Test accuracy force: $force_threshold"
echo "Test accuracy stress: $stress_threshold"
echo "Check accuaracy: $ca"
echo "Test cases file: $cases_file"
echo "Test cases regex: $case"
echo "Generate reference: $g"
echo "--------------------------------"
echo ""


#----------------------------------------------------------
# check_deviation()
#----------------------------------------------------------
check_deviation_pass(){
    deviation=$1
    thr=$2
    echo $(awk -v deviation="$deviation" -v thr="$thr" 'BEGIN{ if (sqrt(deviation*deviation) < thr) print 1; else print 0}')
}

#----------------------------------------------------------
# define a function named 'check_out'
#----------------------------------------------------------
check_out(){
    #------------------------------------------------------
    # input file $1 is 'result.out' in each test directory
    #------------------------------------------------------
    outfile=$1
    thr=$2
    force_thr=$3
    stress_thr=$4
    fatal_thr=$5

    #------------------------------------------------------
    # outfile = result.out
    #------------------------------------------------------
    properties=`awk '{print $1}' $outfile`

    #------------------------------------------------------
    # jd = job description
    #------------------------------------------------------
    if test -e "jd"; then
        jd=`cat jd`
         echo "[----------] $jd"
    fi

    #------------------------------------------------------
    # check every 'key' word
    #------------------------------------------------------
    for key in $properties; do
    
        if [ $key == "totaltimeref" ]; then
            # echo "time=$cal ref=$ref"
            break
        fi

        #--------------------------------------------------
        # calculated value
        #--------------------------------------------------
        cal=`grep "$key" result.out | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # reference value
        #--------------------------------------------------
        ref=`grep "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # computed the deviation between the calculated
        # and reference value
        #--------------------------------------------------
        deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`


        #--------------------------------------------------
        # If deviation < thr, then the test passes,
        # otherwise, the test prints out warning
        # Daye Zheng found bug on 2021-06-20,
        # deviation should be positively defined
        #--------------------------------------------------
        if [ ! -n "$deviation" ]; then
            echo -e "\e[0;31m[ERROR     ] Fatal Error: key $key not found in output.\e[0m"
            let fatal++
            fatal_case_list+=$dir'\n'
            break
        else
            if [ $(check_deviation_pass $deviation $thr) = 0 ]; then
                if [ $key == "totalforceref" ]; then
                    if [ $(check_deviation_pass $deviation $force_thr) = 0 ]; then
                        echo -e "[WARNING   ] "\
                            "$key cal=$cal ref=$ref deviation=$deviation"
                        let failed++
                        failed_case_list+=$dir'\n'
                    else
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                elif [ $key == "totalstressref" ]; then
                    if [ $(check_deviation_pass $deviation $stress_thr) = 0 ]; then
                        echo -e "[WARNING   ] "\
                            "$key cal=$cal ref=$ref deviation=$deviation"
                        let failed++
                        failed_case_list+=$dir'\n'
                    else
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                else
                    echo -e "[WARNING   ] "\
                        "$key cal=$cal ref=$ref deviation=$deviation"
                    let failed++
                    failed_case_list+=$dir'\n'
                fi
                if [ $(check_deviation_pass $deviation $fatal_thr) = 0 ]; then
                    let fatal++
                    fatal_case_list+=$dir
                    echo -e "\e[0;31m[ERROR      ] \e[0m"\
                        "An unacceptable deviation occurs."
                    calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
                    running_path=`echo "OUT.autotest/running_$calculation"".log"`
                    cat $running_path
                fi
                break
            else
                echo -e "\e[0;32m[      OK  ] \e[0m $key"
            fi
        fi
        let ok++
    done
}

#---------------------------------------------
# function to read the threshold from the file
#---------------------------------------------
get_threshold()
{
    threshold_f=$1
    threshold_name=$2
    default_value=$3
    if [ -e $threshold_f ]; then 
        threshold_value=$(awk -v tn="$threshold_name" '$1==tn {print $2}' "$threshold_f")
         if [ -n "$threshold_value" ]; then
            echo $threshold_value
        else
            echo $default_value
        fi
    else
        echo $default_value
    fi
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e $cases_file || (echo "Please specify test cases file by -f option." && exit 1)
which $abacus > /dev/null || (echo "No ABACUS executable was found." && exit 1)

testdir=`cat $cases_file | grep -E $case`
failed=0
failed_case_list=()
ok=0
fatal=0
fatal_case_list=()
fatal_threshold=1
report=""
repo="$(realpath ..)/"

if [ "$sanitize" == true ]; then
    echo "Testing with Address Sanitizer..."
    mkdir ../html
    echo -e "# Address Sanitizer Diagnostics\n" > ../html/README.md
    report=$(realpath ../html/README.md)
    export ASAN_OPTIONS="log_path=asan"
fi

for dir in $testdir; do
    if [ ! -d $dir ];then
        echo -e "\e[0;31m[ERROR     ]\e[0m $dir is not a directory.\n"
        let fatal++
        fatal_case_list+=$dir'\n'
        continue
    fi
    cd $dir
    echo -e "\e[0;32m[ RUN      ]\e[0m $dir"
    TIMEFORMAT='[----------] Time elapsed: %R seconds'
    #parallel test
    time {
        if [ "$case" = "282_NO_RPA" -o "$dir" = "102_PW_BPCG" ]; then
            mpirun -np 1 $abacus > log.txt
        else
            mpirun -np $np $abacus > log.txt
        fi

        if [ "$sanitize" == true ]; then
            echo -e "## Test case ${dir}\n" >> ${report}
            for diagnostic in asan.*; do
                echo -e "### On process id ${diagnostic}\n" >> ${report}
                echo -e "\`\`\`bash" >> ${report}
                cat ${diagnostic} >> ${report}
                echo -e "\`\`\`\n" >> ${report}
            done
        fi
        #$abacus > log.txt
        test -d OUT.autotest || (echo "No 'OUT.autotest' dir presented. Some errors may happened in ABACUS." && exit 1)
        if test -z $g
        then
            ../tools/catch_properties.sh result.out
            if [ $? == 1 ]; then
                echo -e "\e[0;31m [ERROR     ]  Fatal Error in catch_properties.sh \e[0m"
                let fatal++
                fatal_case_list+=$dir'\n'
                break
            else
                my_threshold=$(get_threshold $threshold_file "threshold" $threshold)
                my_force_threshold=$(get_threshold $threshold_file "force_threshold" $force_threshold)
                my_stress_threshold=$(get_threshold $threshold_file "stress_threshold" $stress_threshold)
                my_fatal_threshold=$(get_threshold $threshold_file "fatal_threshold" $fatal_threshold)
                check_out result.out $my_threshold $my_force_threshold $my_stress_threshold $my_fatal_threshold
            fi
        else
            ../tools/catch_properties.sh result.ref
        fi
    }
    echo ""
    cd ../
done

if [ "$sanitize" == true ]; then
    if [[ `uname` == "Darwin" ]]; then
        sed -i '' "s,${repo},,g" ${report}
    else
        sed -i "s,${repo},,g" ${report}
    fi
fi

if [ -z $g ]
then
if [ $failed -eq 0 ]
then
    echo -e "\e[0;32m[ PASSED   ] \e[0m $ok test cases passed."
else
    echo -e "[WARNING]\e[0m    $failed test cases out of $[ $failed + $ok ] failed."
    echo -e $failed_case_list
    if [ $fatal -gt 0 ]
    then
        echo -e "\e[0;31m[ERROR     ]\e[0m $fatal test cases out of $[ $failed + $ok ] produced fatal error."
        echo -e $fatal_case_list
        exit 1
    fi
fi
else
echo "Generate test cases complete."
fi
