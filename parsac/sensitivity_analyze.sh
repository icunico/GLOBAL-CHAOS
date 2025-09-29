#! /bin/bash

pickle_file=fussmann.pickle
pickle_file_analyze=fussmann.analyze.cv.pickle
echo picklefile $pickle_file
file_txt=fussmann_CV.txt
parsac sensitivity analyze $pickle_file --pickle=$pickle_file_analyze cv  > $file_txt
