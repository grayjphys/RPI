#!/bin/bash
for file in * ; do mv $file ${file//MINUS/PLUS} ; done

