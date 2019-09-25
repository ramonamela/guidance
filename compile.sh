#!/bin/bash -e

  # Compile java sources
  mvn clean package

  # Run local Sonar analysis
  # $SONAR_HOME/bin/sonar.sh start
  # mvn sonar:sonar -Dsonar.host.url=http://localhost:9000 -Dsonar.login=44f6a2c96fc322d440e2227bf62d155aa92d8dae

