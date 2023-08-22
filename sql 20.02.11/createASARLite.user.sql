CREATE USER "asar" WITH PASSWORD 'asar' NAME 'ASAR Explorer' SCHEMA "sys";
CREATE SCHEMA "asarlite" AUTHORIZATION "asar";
ALTER USER "asar" SET SCHEMA "asarlite";
CREATE USER "asarpublic" WITH PASSWORD 'asar' NAME 'ASAR Viewer' SCHEMA "sys";
ALTER USER "asarpublic" SET SCHEMA "asarlite";