# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:52:04 2020

@author: admin
"""
# Import biblioteki Biopython oraz zestawu metod PDB.
import Bio
from Bio import PDB

d = PDB.PDBList()
r = PDB.PDBParser()

# Przygotowanie inputa do wprowadzenia kodu PDB.
kod = input("Podaj kod PDB:")

# 1. Pobranie danych z bazy danych na podstawie kody PDB i zapis do pliku.
plik = d.retrieve_pdb_file(kod, pdir='.',file_format='pdb' )

# 2. Ekstrakcja struktury białka z zapisu pliku.
struktura = r.get_structure(kod, plik)

# 3. Wyciągnięcie potrzebnych danych z uzyskanej zmiennej "struktura": 
# nazwa pliku PDB,
# nazwa cząsteczki STR,
# metoda uzyskania struktury,
# rozdzielczosc struktury białka,
# czasopismo publikacji,
# brakujące aminokwasy.

nazwaPDB      = struktura.header['name'].title()
nazwaSTR      = struktura.header['compound']['1']['molecule']  
metoda        = struktura.header['structure_method']
rozdzielczosc = struktura.header['resolution']
czasopismo    = struktura.header['journal_reference']
missing_residues  = struktura.header['missing_residues']

# 4. Iteracja po obiekcie "struktura" i uzyskanie l.rezuidów, l.cząsteczek wody, l.ligandów oraz wyliczenie odległosci pomiędzy pierwszym i ostatnim CA.
Lamino = 0
Lwody  = 0
Llig   = 0
a1 = struktura[0]['A'].child_list[0]['CA']
for r in struktura.get_residues():
    if r.id[0] == ' ':
        Lamino +=1
        a2 = r['CA']
    elif r.id[0] == 'W':
        Lwody += 1
    else:
        Llig +=1  

# 4. Wyliczenie brakujących aminokwasów.
listaBrakujacychAminokwasow = []
for r in missing_residues:
    listaBrakujacychAminokwasow.append(r['res_name'])
        
# Zwrot uzyskanych rezultatów: 
    
print("\nPlik o kodzie '%s' to:\n\t%s" % (kod,nazwaPDB.upper()))
print("\nPlik ten zawiera czasteczke o nazwie: '%s'" % (nazwaSTR))
print ("\nMetoda uzyskania struktury to: %s" % (metoda.upper()))

if "X-RAY" in metoda.upper():
    print ("\nRozdzielczos struktury: %s" % (rozdzielczosc))

elif "NMR" in metoda.upper():
    print ("\nStruktura zawiera %s modeli" % (len(struktura)))
    
else:
    pass

print ("\nOpublikowana zostala w: \n\t%s" %(czasopismo))

print("\nLiczba reziduów w strukturze: %d" % (Lamino))
print("\nLiczba czasteczek wody      : %d" % (Lwody))
print("\nLiczba ligandow w strukturze: %d" % (Llig))
print('\nLiczba brakujących aminokwasów: %s' % (len(listaBrakujacychAminokwasow)))
print('Są to kolejno: %s' %  (listaBrakujacychAminokwasow))        
print("\nOdleglosc pierwszego i ostatniego CA: %s\n\n\n" % (a1-a2))