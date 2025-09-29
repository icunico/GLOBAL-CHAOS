#! /bin/bash

pickle_file=fussmann.pickle
xml_file=fussmann.xml
echo xmlfile $xml_file
echo picklefile $pickle_file
parsac sensitivity sample $xml_file $pickle_file random 100000
