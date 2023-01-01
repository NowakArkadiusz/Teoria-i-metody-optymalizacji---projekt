from tkinter import *
import projektv10_funkcje as funkcje

root = Tk()
root.title("Aplikacja demo- jeszcze :p")


def pobierz_dane_z_interfejsu():

    maksymalny_numer = 0
    p=0
    j=0
    l=0
    punkt_poczatkowy = []
    znak_niewiadomej = []
    sczytana_funkcja = []
    F = pole_funkcja.get()
    punkt_poczatkowy_sczytany = pole_punkt_poczatkowy.get()

    # Sczytanie punktow poczatkowych do tablicy
    
    for i in range(int(len(punkt_poczatkowy_sczytany))):
        if ( punkt_poczatkowy_sczytany[i] == ';' ) or ( punkt_poczatkowy_sczytany[i] == ',' ):
            punkt_poczatkowy.append(punkt_poczatkowy_sczytany[i-j:i])
            j=0
        else:
            j=j+1

       
    # ograniczenia dla calego algorytmu
    epsilon = float(pole_eps.get())
    epsilon1 = float(pole_eps1.get())
    epsilon2 = float(pole_eps2.get())
    L = int(pole_l.get())

    # ograniczenia dla metody zlotego podzialu
   
    #m = float(pole_m_metoda.get())
    #epsilon_metody = float(pole_eps_metoda.get())
    #L_zlotego_podzialu = int(pole_l_metoda.get())


    # ograniczenia dla metody bisekcji z testem dwuskosnym glodensteina

    poczatkowe_tau_r =float(pole_tau_r.get())
    beta = float(pole_beta.get())
    
   
   # Sprawdzanie ilosci niewiadomoych i symboli niewiadomych
   # 1) Jezeli dany znak to litera i tylko jedna litera (bo inaczej to moze byc funkcja typu sin, cos, pi, exp, itd) to:
   # 2) Spradz czy to pierwsza niewiadoma, Jezeli tak to podstaw do tablicy niewiadomych i zlicz pierwsa niewiadoma
   # 3) Nastepnie sprawdz czy dana niewiadoma juz nie znajduje sie w tablicy niewiadomych, jezeli sie znajduje przerwij petle i sprawdz kolejny znak (wroc do punktu nr 1)
   # 4) Jezeli dana litera nie znajduje sie w tablicy niewiadomych podstaw ja do tejze tablicy i zlicz kolejna niewiadoma do zmiennej p
   # 5) Wroc do punktu nr 5 i powtarzaj tak az nie sprawdisz wszystkich wprowadzonych znakow


    for i in range(int(len(F))-1):    
      
        if F[i].isalpha()==True and F[i-1].isalpha()==False and F[i+1].isalpha()==False:            #1)
            if l==0:                                                                                #2)
                znak_niewiadomej.append(F[i])
                l=l+1
                p=p+1
            for k in range(l):
                if znak_niewiadomej[k]==F[i]:                                                       #3)
                    break
                elif k==l-1:                                                                        #4)
                    znak_niewiadomej.append(F[i])
                    l=l+1
                    p=p+1
        if i==int(len(F)-2) and F[i].isalpha()==False and F[i+1].isalpha()==True:                   # sprawdzenie ostatniego znaku
            if l==0:                                                                                #2)
                znak_niewiadomej.append(F[i+1])
                l=l+1
                p=p+1
            for k in range(l):
                if znak_niewiadomej[k]==F[i+1]:                                                     #3)
                    break
                elif k==l-1:                                                                        #4)
                    znak_niewiadomej.append(F[i+1])
                    l=l+1
                    p=p+1
        if i==0 and F[i].isalpha()==True and F[i+1].isalpha()==False:                               # sprawdzenie pierwszego znaku
            if l==0:                                                                                #2)
                znak_niewiadomej.append(F[i])
                l=l+1
                p=p+1
            for k in range(l):
                if znak_niewiadomej[k]==F[i]:                                                     #3)
                    break
                elif k==l-1:                                                                        #4)
                    znak_niewiadomej.append(F[i])
                    l=l+1
                    p=p+1
            
   
    #return  punkt_poczatkowy, epsilon, epsilon1, epsilon2, epsilon_metody, m, L_zlotego_podzialu, L, p, F,znak_niewiadomej

    return  punkt_poczatkowy, epsilon, epsilon1, epsilon2, poczatkowe_tau_r, beta, L, p, F,znak_niewiadomej


def przyblizenie(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p):
    wykres_od = float(pole_wykres_od.get())
    wykres_do = float(pole_wykres_do.get())
    dokladnosc = 1000

    wykres_od = wykres_od/2
    wykres_do = wykres_do/2

    pole_wykres_od.delete(0,END)
    pole_wykres_do.delete(0,END)
    pole_wykres_od.insert(END,wykres_od)
    pole_wykres_do.insert(END,wykres_do)

    funkcje.rysowanie_wykresu(wykres_od, wykres_do,1000,F,X,wartosc_funkcji, k, p) 

def oddalenie (pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p):
    wykres_od = float(pole_wykres_od.get())
    wykres_do = float(pole_wykres_do.get())
    dokladnosc = 1000

    wykres_od = wykres_od*2
    wykres_do = wykres_do*2

    pole_wykres_od.delete(0,END)
    pole_wykres_do.delete(0,END)
    pole_wykres_od.insert(END,wykres_od)
    pole_wykres_do.insert(END,wykres_do)

    funkcje.rysowanie_wykresu(wykres_od, wykres_do,1000,F,X,wartosc_funkcji, k, p) 

def centruj (pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p):
    wykres_od = abs(float(pole_wykres_od.get()))
    wykres_do = abs(float(pole_wykres_do.get()))
    dokladnosc = 1000

    if wykres_od>= wykres_do:
        wykres_do=wykres_od
        wykres_od=-1*wykres_od     
    else:
        wykres_od=-1*wykres_do
        wykres_do=wykres_do  

    pole_wykres_od.delete(0,END)
    pole_wykres_do.delete(0,END)
    pole_wykres_od.insert(END,wykres_od)
    pole_wykres_do.insert(END,wykres_do)

    funkcje.rysowanie_wykresu(wykres_od, wykres_do,1000,F,X,wartosc_funkcji, k, p) 

def wykres(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p):
    wykres_od = float(pole_wykres_od.get())
    wykres_do = float(pole_wykres_do.get())
    dokladnosc = 1000
    funkcje.rysowanie_wykresu(wykres_od, wykres_do,1000,F,X,wartosc_funkcji, k, p) 



def oblicz():
 
    #punkt_poczatkowy, epsilon, epsilon1, epsilon2, epsilon_metody, m, L_zlotego_podzialu, L, p, F, znak_niewiadomej = pobierz_dane_z_interfejsu()
    #X, wartosc_funkcji, k, kryterium_stopu, wartosc_kryterium_stopu = funkcje.metoda_polaka_ribierry (punkt_poczatkowy, epsilon, epsilon1, epsilon2, epsilon_metody, m, L_zlotego_podzialu, L, p, F, znak_niewiadomej)

    punkt_poczatkowy, epsilon, epsilon1, epsilon2, poczatkowe_tau_r, beta_gold, L, p, F, znak_niewiadomej = pobierz_dane_z_interfejsu()
                                         
    X, wartosc_funkcji, k, kryterium_stopu, wartosc_kryterium_stopu = funkcje.metoda_polaka_ribierry (punkt_poczatkowy, epsilon, epsilon1, epsilon2, poczatkowe_tau_r, beta_gold, L, p, F, znak_niewiadomej)

    # przyblizenie wartosci 
    for i in range(k+1):
        wartosc_funkcji[i] = round(wartosc_funkcji[i],4)
        for j in range(p):
            X[i][j] = round(X[i][j],4)


    root2 = Toplevel(root)
    root2.title("Drugie okieneczko")
    

    # definiowanie obramowań   
    obramowanie_rozwiazanie = LabelFrame(root2, text="Rozwiazanie: ") 
    obramowanie_kolejne_iteracje = LabelFrame(root2, text="Kolejne iteracje: ")
    obramowanie_wykres = LabelFrame(root2, text="Wykres ")

    # polozenie obramowan na ekranie
    obramowanie_rozwiazanie.grid(row=0,column=0)
    obramowanie_kolejne_iteracje.grid(row=1,column=0)
    obramowanie_wykres.grid(row=0,column=1)

    # definicja wyświetlonych tekstów
    wsp_pkt = ""
    pkt = ""
    for i in range(int(len(znak_niewiadomej))):
        pkt = pkt + str(znak_niewiadomej[i]) + ","
        wsp_pkt = wsp_pkt + str(znak_niewiadomej[i]) + " = " + str(X[k][i]) + "\n"
    pkt=pkt[0:int(len(pkt))-1]

    tekst_wartosc_funkcji = Label(obramowanie_rozwiazanie ,text = "Wartość funkcji: \n" +"f(" + str(pkt) + ") = " + str(wartosc_funkcji[k])) 
    tekst_wspolrzedne_punktu = Label(obramowanie_rozwiazanie ,text = "Współrzędne punktu: \n"+ wsp_pkt)  
    tekst_kryterium_stopu = Label(obramowanie_rozwiazanie ,text = "Kryterium stopu:     " + str(kryterium_stopu) + " \n")
    tekst_wartosc_kryterium_stopu = Label(obramowanie_rozwiazanie ,text = "Wartość kryterium stopu:     " + str(wartosc_kryterium_stopu) + " \n")

    # definicja tablicy oraz wyswietlanie na niej kolejnych iteracji
    label_tekst = Text(obramowanie_kolejne_iteracje, width =60, height = 40, wrap = WORD)     # WRAP = WORD - cos jak enter jak zabraknie linijki
    for i in range(k+1):
        tekst = str('Krok nr ' + str(i) + ': f(' + str(X[i]) + ') = ' + str(wartosc_funkcji[i]) + '\n')
        label_tekst.insert(END, tekst)


    #definicja pola wejsciowego

    pole_wykres_od = Entry(obramowanie_wykres, width=10)
    pole_wykres_do = Entry(obramowanie_wykres, width=10)

    if p==1:
        if X[k][0] >= wartosc_funkcji[k]:
            pole_wykres_od.insert(END, wartosc_funkcji[k] - abs( (X[k][0]+wartosc_funkcji[k])/2 ))
            pole_wykres_do.insert(END, X[k][0] + abs( (X[k][0]+wartosc_funkcji[k])/2 ))
        else:
            pole_wykres_od.insert(END, X[k][0] - abs( (X[k][0]+wartosc_funkcji[k])/2 ))
            pole_wykres_do.insert(END, wartosc_funkcji[k] + abs( (X[k][0]+wartosc_funkcji[k])/2 ))

    if p==2:
        if X[k][0] >= X[k][1]:
            pole_wykres_od.insert(END, X[k][1] - abs( (X[k][0]+X[k][1])/2 ))
            pole_wykres_do.insert(END, X[k][0] + abs( (X[k][0]+X[k][1])/2 ))
        else:
            pole_wykres_od.insert(END, X[k][0] - abs( (X[k][0]+X[k][1])/2 ))
            pole_wykres_do.insert(END, X[k][1] + abs( (X[k][0]+X[k][1])/2 ))

    # definicja przyciskow

    przycisk_oddalenie = Button(obramowanie_wykres, text = "Oddalenie", command= lambda: oddalenie(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p))
    przycisk_przyblizenie = Button(obramowanie_wykres, text = "Przybliżenie", command= lambda: przyblizenie(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p))
    przycisk_wykres = Button(obramowanie_wykres, text = "Rysuj wykres", command= lambda: wykres(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p))
    przycisk_centruj = Button(obramowanie_wykres, text = "Centruj", command= lambda: centruj(pole_wykres_od,pole_wykres_do,F,X,wartosc_funkcji, k, p))

    #polozenie elementow na ekranie 
    tekst_wartosc_funkcji.grid(row=0,column=0)
    tekst_wspolrzedne_punktu.grid(row=1,column=0)
    tekst_kryterium_stopu.grid(row=2,column=0)
    tekst_wartosc_kryterium_stopu.grid(row=3,column=0)

    label_tekst.grid(row=0,column=0)

    przycisk_oddalenie.grid(row=0,column=0)
    przycisk_przyblizenie.grid(row=0,column=1)
    przycisk_wykres.grid(row=0,column=2)
    pole_wykres_od.grid(row=1,column=0)
    pole_wykres_do.grid(row=1,column=1)
    przycisk_centruj.grid(row=1,column=2)


def funkcja1():
    pole_funkcja.delete(0,END)
    pole_funkcja.insert(END,"(x**2 + y - 11)**2 + (x + y**2 -7 )**2 - 200 ")  #(1.5-x+x*y)**2+(2.25-x+x*y**2)**2+(2.625-x+x*y**3)**2 punkt 1;1; do niesk.(gradient wynosi zero) zmien druga niewiadoma na np.1.5      #(x+2*y-7)**2+(2*x+y-5)**2  # x**2+x*y+0.5*y**2-x-y
    pole_punkt_poczatkowy.delete(0,END)                                             #x**2+x*y-0.5*y**2-x-y   #sin(x+y)+(x-y)**2-1.5*x+2.5*y+1       # -20*exp(-0.2*sqrt(0.5*(x**2+y**2)))-exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+exp(1)+20
    pole_punkt_poczatkowy.insert(END,"1;1;")                                       # 0.26*(x**2+y**2)-0.48*x*y    #(sin(3*pi*x))**2+(x-1)**2*(1+(sin(3*pi*y))**2)+(y-1)**2*(1+(sin(2*pi*y))**2)
    pole_eps.delete(0,END)                                                         # 6 - hump camel (4-2.1*x**2+x**4/3) * x**2 + x*y + (-4 + 4 *y**2) * y**2
    pole_eps.insert(END,"0.001")                                                   # 3 - hump camel 2*x**2-1.05*x**4+x**6/6 + x*y+y**2
    pole_eps1.delete(0,END)                                                        # x**4+y**4-x**2-y**2               # x**2+2*y**2-6*x+x*y usun pozniej, to z wykladu
    pole_eps1.insert(END,"0.001")                                                  # (x-2)**2+(y-2)**2    punkt pooczatkowy 1;1;     kolejny punkt powininen byc od razu 2;2;        # (x-4)**2 + (x-y**2)**2 funckja 
    pole_eps2.delete(0,END)
    pole_eps2.insert(END,"0.001")
    pole_l.delete(0,END)
    pole_l.insert(END,"30")
    #pole_m_metoda.delete(0,END)
    #pole_m_metoda.insert(END,"1") # 0.001
    #pole_eps_metoda.delete(0,END)
    #pole_eps_metoda.insert(END,"0.0001")
    #pole_l_metoda.delete(0,END)
    #pole_l_metoda.insert(END,"1000")
    pole_tau_r.delete(0,END)
    pole_tau_r.insert(END,"9")
    pole_beta.delete(0,END)
    pole_beta.insert(END,"0.4")

def funkcja2():
    pole_funkcja.delete(0,END)
    pole_funkcja.insert(END,"4*x**2 - 2.1 * x**4 + (1/3) * x**6 + x*y -4 * y**4")

def funkcja3():
    pole_funkcja.delete(0,END)
    pole_funkcja.insert(END,"x**2 -x +1")

def clear():
    pole_funkcja.delete(0,END)
    pole_punkt_poczatkowy.delete(0,END)
    pole_eps.delete(0,END)
    pole_eps1.delete(0,END)
    pole_eps2.delete(0,END)
    pole_l.delete(0,END)
    #pole_m_metoda.delete(0,END)
    #pole_eps_metoda.delete(0,END)
    #pole_l_metoda.delete(0,END)
    pole_tau_r.delete(0,END)
    pole_beta.delete(0,END)

# definiowanie obramowań   
obramowanie_funkcji = LabelFrame(root, text="Funkcja", padx=10, pady=10)
obramowanie_punkt_poczatkowy = LabelFrame(root, text="Punkt poczatkowy", padx=30, pady=10)
obramowanie_kryteria_stopu_algorytmu = LabelFrame(root, text="Kryteria stopu algorytmu", padx=114)
obramowanie_kryteria_stopu_metody = LabelFrame(root, text="Kryteria stopu metody bisekcji z testem dwus. Gold.", padx=100, pady=10)
obramowanie_przykladowe_funkcje = LabelFrame(root, text="Przykladowe funkcje", padx=50, pady=20)
obramowanie_przyciski = LabelFrame(root, text="Przyciski funkcyjne", padx=40, pady=6)     # JAK INACZEJ TO NAZWAC ????????????????



# polozenie obramowan na ekranie

obramowanie_funkcji.grid(row=0,column=0, padx=5, pady=5)
obramowanie_punkt_poczatkowy.grid(row=0,column=1, padx=5, pady=5)
obramowanie_kryteria_stopu_algorytmu.grid(row=1,column=0, padx=5, pady=5)
obramowanie_kryteria_stopu_metody.grid(row=1,column=1, padx=5, pady=5)
obramowanie_przykladowe_funkcje.grid(row=2,column=0, padx=5, pady=5) 
obramowanie_przyciski.grid(row=2,column=1, padx=5, pady=5) 


# definiowanie pol do wpisywania funkcji/liczba/danych przez uzytkownika
pole_funkcja = Entry(obramowanie_funkcji, width=50) # , width=50

pole_punkt_poczatkowy = Entry(obramowanie_punkt_poczatkowy, width=38)

pole_eps = Entry(obramowanie_kryteria_stopu_algorytmu, width=10)
pole_eps1 = Entry(obramowanie_kryteria_stopu_algorytmu, width=10)
pole_eps2 = Entry(obramowanie_kryteria_stopu_algorytmu, width=10)
pole_l =  Entry(obramowanie_kryteria_stopu_algorytmu, width=10)

#pole_m_metoda = Entry(obramowanie_kryteria_stopu_metody, width=10)
#pole_eps_metoda = Entry(obramowanie_kryteria_stopu_metody, width=10)
#pole_l_metoda = Entry(obramowanie_kryteria_stopu_metody, width=10)

pole_tau_r = Entry(obramowanie_kryteria_stopu_metody, width=10)
pole_beta = Entry(obramowanie_kryteria_stopu_metody, width=10)




# definicja wyświetlonych tekstów

tekst_eps = Label(obramowanie_kryteria_stopu_algorytmu ,text = "Eps:")
tekst_eps1 = Label(obramowanie_kryteria_stopu_algorytmu ,text = "Eps1:")
tekst_eps2 = Label(obramowanie_kryteria_stopu_algorytmu ,text = "Eps2:")
tekst_l = Label(obramowanie_kryteria_stopu_algorytmu ,text = "L:")

#tekst_m_metoda = Label(obramowanie_kryteria_stopu_metody ,text = "M:")
#tekst_eps_metoda = Label(obramowanie_kryteria_stopu_metody ,text = "Eps:")
#tekst_l_metoda = Label(obramowanie_kryteria_stopu_metody ,text = "L:")

tekst_tau_r = Label(obramowanie_kryteria_stopu_metody ,text = "Tau_R:")
tekst_beta = Label(obramowanie_kryteria_stopu_metody ,text = "Beta:")

# definiowanie przyciskow

przycisk_funkcja1 = Button(obramowanie_przykladowe_funkcje, text = "Funkcja nr.1", command=funkcja1) # padx=40, pady=10,
przycisk_funkcja2 = Button(obramowanie_przykladowe_funkcje, text = "Funkcja nr.2", command=funkcja2)
przycisk_funkcja3 = Button(obramowanie_przykladowe_funkcje, text = "Funkcja nr.3", command=funkcja3)

przycisk_clear = Button(obramowanie_przyciski, text = "Wyczyść formularz",  command=clear, padx=10, pady=10) # padx=10, pady=10,
przycisk_oblicz = Button(obramowanie_przyciski, text = "Oblicz",  command=oblicz, padx=10, pady=10) # padx=10, pady=10,

#polozenie elementow na ekranie 

pole_funkcja.grid(row=0,column=0)

pole_punkt_poczatkowy.grid(row=0,column=0)

tekst_eps.grid(row=0,column=0)
pole_eps.grid(row=0,column=1)
tekst_eps1.grid(row=1,column=0)
pole_eps1.grid(row=1,column=1)
tekst_eps2.grid(row=2,column=0)
pole_eps2.grid(row=2,column=1)
tekst_l.grid(row=3,column=0)
pole_l.grid(row=3,column=1)

#tekst_m_metoda.grid(row=0,column=0)
#pole_m_metoda.grid(row=0,column=1)
#tekst_eps_metoda.grid(row=1,column=0)
#pole_eps_metoda.grid(row=1,column=1) 
#tekst_l_metoda.grid(row=2,column=0)
#pole_l_metoda.grid(row=2,column=1)

tekst_tau_r.grid(row=0,column=0)
pole_tau_r.grid(row=0,column=1)
tekst_beta.grid(row=1,column=0)
pole_beta.grid(row=1,column=1)


przycisk_funkcja1.grid(row=0,column=0)
przycisk_funkcja2.grid(row=0,column=1)
przycisk_funkcja3.grid(row=0,column=2)

przycisk_clear.grid(row=0,column=0,padx=5, pady=5)
przycisk_oblicz.grid(row=0,column=1,padx=5, pady=5)


root.mainloop()     # loopuje cala funkcje (jest po to gdyz program sprawdza polozenie kursora wykonuje petle, sprawdza polozenie itd. wkolko, przerwanie tej petli konczy dzialanie programu)



