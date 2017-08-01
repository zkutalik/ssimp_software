# This project, ssimp, doesn't really use 'make'. It was developed with 'stu',
# an excellent replacement for 'make'. See https://github.com/kunegis/stu for
# more. In this simple Makefile, we create the stu binary and then use it to
# build the ssimp binary in 'bin/ssimp'

bin/ssimp: Makefile stu/stu
	stu/stu bin/ssimp

stu/stu:
	set -e; cd stu; ls
	set -e; cd stu; ${CXX} stu.cc -o stu -std=c++11
