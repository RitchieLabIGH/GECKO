g++ -std=c++0x -O3 -fopenmp -o indexingKmers IndexingKMers.cpp
g++ -std=c++0x -O3 -fopenmp RemoveNonInformativeKmers.cpp -o NIK
g++ -std=c++0x -O3 -fopenmp JoinMultiImportationFiles.cpp -o JoinKMers
g++ -std=c++0x -O3 -fopenmp JoinUltraBigMatrix.cpp -o JoinKmersUltra
g++ -std=c++0x -O3 -fopenmp transformIntoMLformat.cpp -o transformIntoML
g++ -std=c++0x -O3 -fopenmp CutMatrixByClass.cpp -o CutMatrixByClass
g++ -std=c++0x -O3 MLformatToBinary.cpp -lm -lz -o transformIntoBinary
g++ -std=c++0x -O3 MLformatToBinaryRemoveGroup.cpp -lm -lz -o transformIntoBinaryRemoveGroup
make -f ./anova_makefile
mv anovaFilter indexingKmers NIK JoinKmersUltra JoinKMers transformIntoML transformIntoBinary transformIntoBinaryRemoveGroup CutMatrixByClass ../.
