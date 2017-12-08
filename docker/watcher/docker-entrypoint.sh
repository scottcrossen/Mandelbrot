#!/bin/bash

SLEEPAMNT=2

function finish {
	if [ $PID1 -ge 0 ]; then
		kill $PID1
	fi
	kill $PID2
}

trap 'finish' SIGTERM

function test {
	echo -e "\nTesting Started."
	if [ ! -d "./test" ]; then
		mkdir test
	fi
	./rel/mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 8192 test/pic.ppm
	convert -normalize test/pic.ppm test/pic.png
	echo -e "\nTesting Finished"
}

if [ ! -f ./Makefile ]; then
	echo -e "\n[INFO] No Makefile found. File must exist\n"
	exit
fi

PID1=-1
while true; do
	if [ $PID1 -ge 0 ]; then
		echo -e "\n[INFO] Restarting test cases in ${SLEEPAMNT} seconds\n"
		kill $PID1
		PID1=-1
	fi
	sleep ${SLEEPAMNT} && test &
	PID1=$!
	inotifywait -e modify -e delete -e create -e attrib --exclude ./test ./*
	make
done &
PID2=$!

wait
