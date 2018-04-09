### AUTHOR: Spencer Richman ###
### CONTACT: scr6045@rit.edu ###

# This script serves as an overly complicated wrapped for various Entrez eutils
# Esearch for biosample Ids based on criterion specified in file --> Epost in batches of 1000 ids -->
# --> elink each batch for linked SRA archives --> return links --> archive into tar.gz files

## THIS IS NOT USER-FRIENDLY NOR IS IT EASILY EXTENSIBLE IN ITS CURRENT FORM 
