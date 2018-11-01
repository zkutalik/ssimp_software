# This project, ssimp, doesn't really use 'make'. It was developed with 'stu',
# an excellent replacement for 'make'. See https://github.com/kunegis/stu for
# more. In this simple Makefile, we create the stu binary and then use it to
# build the ssimp binary in 'bin/ssimp'

.PHONY: stu/stu bin/ssimp

bin/ssimp: Makefile stu/stu
	@echo
	@echo "   Building 'bin/ssimp'"
	@echo
	stu/stu bin/ssimp

stu/stu:
	@echo
	@echo "   Compiling 'stu', which is then used to build 'ssimp'"
	@echo "   Requires g++ , the GNU C++ compiler"
	@echo
	set -e; cd stu; ${CXX} stu.cc -o stu -std=c++11
