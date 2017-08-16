
clone ssimp_software
cp ssimp_software ssimp_linux_{x}
cd ssimp_linux_{x}
rm .git file
rm -r compiled/
stu bin/ssimp
cd ..
zip -r ssimp_linux_x{.zip,}
