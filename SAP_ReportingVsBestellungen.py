# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 09:48:48 2022

@author: just_d
"""
import pandas as pd
import numpy as np

"""
Configuration
Angabe zum Dateipfad
pwdBest Dateipfad zum File aus den Bestellungen
PwdKost Dateipfad zum File aus dem Projektkostenreport
"""
pwdBest = r"C:\Users\just_d\Desktop\tmp\SAP\SAP 8.16004.003.01.01 Athos Consolidaion Athos Beamline Bestellungen.xlsx"
pwdKost = r"C:\Users\just_d\Desktop\tmp\SAP\SAP 8.16004.003.01.01 Athos Consolidaion Athos Beamline Projektkosten.xlsx"

"""
Configuration ends here
"""


best = pd.read_excel(pwdBest)
kost = pd.read_excel(pwdKost)

kost_einkaufsbeleg = kost[["Betrag in BukrsWähr.", "Bestellung"]].groupby("Bestellung").sum()
best_einkaufsbeleg = best[["Betrag in Buchungskreiswährung", "Einkaufsbeleg"]].groupby("Einkaufsbeleg").sum()


kost_einkaufsbeleg = kost_einkaufsbeleg.sort_values("Bestellung")
best_einkaufsbeleg = best_einkaufsbeleg.sort_values("Einkaufsbeleg")

i = 0
s = 0
for i in best["Betrag in Buchungskreiswährung"]:
    s = s+i
print('Summe aus Bestellungen ist: ' + str(round(s,2)) +' CHF')
sumBest = s    
i = 0
j= 0
s = 0

for i in kost["Betrag in BukrsWähr."]:
    s= s+i
print('Summe aus Kosten ist: ' +str(round(s,2)) +' CHF')
sumKost = s 
print('Unterschied Kosten aus Reporting zu Bestellung ' +str(round(sumKost-sumBest,2)) +' CHF')
print()
#Vergleiche Bestellungen mit Reporting
#index der bestellten ist länger als der index dejenigen welcher Kosten verurscahen deshalb wird geprüft wie viele das sind
bestNoCost = []
comp  = []
for ind in best_einkaufsbeleg.index:
    try:
        kost_einkaufsbeleg.loc[[ind]] #shaut ob index Bsp. 450001336 vorhanden ist
        dif = kost_einkaufsbeleg.loc[[ind]].iloc[0,0] - best_einkaufsbeleg.loc[[ind]].iloc[0,0]
        comp.append([ind,round(dif,2)])
    except:
        bestNoCost.append([ind,best_einkaufsbeleg.loc[[ind]].iloc[0,0]])
        #print(best_einkaufsbeleg.loc[[ind]].iloc[0,0])
       
    
comp = pd.DataFrame(comp, columns = ['Best. Nr.','Dif. Bestellt / Kosten CHF'])
comp=comp.sort_values(by=['Dif. Bestellt / Kosten CHF'],ascending=False)
print('Abweichungen Bestellt zu Kosten mit der selben Bestellnummer ' +str(comp.sum()[1]) +'CHF')
print('Die fünf grösten Punkte sind ')
print(comp.head(5))
print()


NoBest=kost[kost['Bestellung'].isnull()] #gibt alle Zeilen aus dern Bestellung Feld leer ist
NoBest=NoBest.sort_values('Betrag in BukrsWähr.',ascending=False)
sumNoBest=NoBest.sum(numeric_only=True)[1]
print('Kosten ohne Bestellung ' +str(round(sumNoBest,2)) + 'CHF')
print('Die zehn grösten Posten sind: ')
print(NoBest[["Betrag in BukrsWähr.","Positionstext","Buchungsdatum", "Belegnummer"]].head(10))

if bestNoCost != []:
    i = 0
    s = 0
    for i in bestNoCost: #berechent die Summe der Bestellten Werte welche nicht im Reporting auftauchen
        s = s + i[1]
    print("Die gesammt Summe Obligo ist " + str(round(s,2)) +' CHF')  
    obligo=s
    print("Folgende Einkaufsbelege wurden Bestellt sind aber nicht in den Kosten, (Obligo)")
    bestNoCost = pd.DataFrame(bestNoCost,columns = ['Best. Nr.','Veranschlagte Kosten CHF']) #bring to pandas format
    bestNoCost =bestNoCost.sort_values('Veranschlagte Kosten CHF',ascending=False)
    print(bestNoCost)
    

tot= sumKost-sumBest
print('Von dem total ' +str(round(tot,2)) +' CHF Unterschied zwischen Reporting und Bestellungen sind ' +str(round(sumNoBest*100/tot,0)) +' % aus Quellen die nicht in Bestellungen auftauchen. ' +str(round(comp.sum()[1]*100/tot,0)) +' % kommen daher, dass mehr abgebucht wurde als in der Bestellung angegeben und ' +str(round(obligo*100/tot*-1,0)) +' % stammen aus Obligo.')

# a =kost_einkaufsbeleg.loc[[4500001336.0]] #gibt Zeile mit index 4500013356 aus
# a.iloc[0,0] #¸gibt wert in spalte null zeile null von a aus
