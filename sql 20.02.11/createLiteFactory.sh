#!/bin/bash

monetdbd create $DB_PATH #/Users/lptolik/Documents/Projects/monetDB/asarLite/
monetdbd start $DB_PATH
monetdb create asar
monetdb release asar
mclient -u monetdb -d asar ./createASARLite.user.sql
mclient -u asar -d asar ./ASARLite.sql

