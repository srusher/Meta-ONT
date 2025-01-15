#!/bin/bash

while IFS= read -r line; do

    if [[ $(echo $line | awk '{print $5}') == "strain" ]]; then

        tax_id=$(echo $line | awk '{print $1}')
        parent_id=$(echo $line | awk '{print $3}')

        sed -i "s/\b$tax_id\b/$parent_id/g" $2

    fi

done < $1