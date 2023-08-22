#!/bin/bash

monetdb stop asar
monetdb lock asar
monetdbd stop $DB_PATH #/Users/lptolik/Documents/Projects/monetDB/asarLite/
rm -rf "$DB_PATH/*"