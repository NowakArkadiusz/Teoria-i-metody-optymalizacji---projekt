import numpy as np
import sympy as sym 
import math
from sympy import lambdify
from sympy import FunctionMatrix, Matrix
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def kierunek_poszukiwan_MN(k, gradient,p,beta,d): 

    for i in range(p):
        d[k][i] = ( - gradient[k][i] )
    print("wartosc d= " + str(d[k]))
    return d

def kierunek_poszukiwan_zle(k, gradient,p,beta,d): 

    if k==0:
        for i in range(p):
            d[k][i] = ( - gradient[k][i] )
    else:      
        for i in range(p):
            gradient_T1_0 = np.transpose(gradient[k][i]-gradient[k-1][i])
            gradient_T2_0 = np.transpose(gradient[k-1][i])    
            beta[k-1][i] = ( gradient_T1_0*gradient[k][i] ) / (gradient_T2_0 * gradient[k-1][i])   
            d[k][i] = ( - gradient[k][i] + beta[k-1][i] * d[k-1][i] )
    
    if k!=0:
        print('wartosc beta = ' + str(beta[k-1]))
    print("wartosc d= " + str(d[k]))

    return d

def kierunek_poszukiwan(k, gradient,p,beta,d): 

    if k==0:
        for i in range(p):
            d[k][i] = ( - gradient[k][i] )
    else:
        gradient_T1_0 = np.transpose(np.subtract(gradient[k],gradient[k-1]))
        gradient_T2_0 = np.transpose(gradient[k-1])    
        print('licznik = ' + str(np.dot(gradient_T1_0,gradient[k])))
        print('mianownik = ' + str(np.dot(gradient_T2_0,gradient[k-1])))
        beta[k-1] = ( np.dot(gradient_T1_0,gradient[k])) / (np.dot(gradient_T2_0,gradient[k-1]))   
        for i in range(p):  
            d[k][i] = ( - gradient[k][i] + beta[k-1] * d[k-1][i] )

    if k!=0:
        print('wartosc beta = ' + str(beta[k-1]))
    print("wartosc d= " + str(d[k]))

    return d
  

def metoda_zlotego_podzialu(epsilon_metody, m, L_zlotego_podzialu, X, d, k, p,f):

    a = [0] * L_zlotego_podzialu
    b = [0] * L_zlotego_podzialu
    tmp = [0] * L_zlotego_podzialu     
    l = [0] * L_zlotego_podzialu       
    v = [0] * L_zlotego_podzialu
    w = [0] * L_zlotego_podzialu
    alfa = [0] * L_zlotego_podzialu
    for i in range(L_zlotego_podzialu):
        a[i] = [0] * p
        b[i] = [0] * p
        tmp[i] = [0] * p
        l[i] = [0] * p
        v[i] = [0] * p
        w[i] = [0] * p
        alfa[i] = [0] * p


    zmienne1= [0] * p
    zmienne2= [0] * p
    zmienne3= [0] * p
    zmienne4= [0] * p              
    zmienne5= [0] * p
    

 
    for i in range(p):
            zmienne1[i] = X[k][i]+ m * d[k][i]

    if f(zmienne1) > f(X[k]):  
        for i in range(p):
            a[0][i] = 0
            b[0][i] = m    
            #print('m jest rowne = ' + str(m))
            #print('pierwszy punkty - ' + str(zmienne1))
    else:
        j = 0 
        for i in range(p):
            tmp[j][i] = 0

        while True:
            if j+1 == L_zlotego_podzialu:                                             
                print("nie znalazlo min. - dazylo caly czas do - niesk.")
                break

            for i in range(p):
                tmp[j+1][i] = tmp[j][i] + m
                zmienne2[i] = X[k][i] + tmp[j+1][i]*d[k][i]
                zmienne3[i] = X[k][i] + tmp[j][i]*d[k][i]
               # print('kolejne punkty - ' + str(zmienne2))

            if f(zmienne2) >= f(zmienne3): 
              
                for i in range(p):
                    a[0][i] = tmp[j-1][i]   
                    b[0][i] = tmp[j+1][i]               
                break
            else:
                j=j+1
                   
    i=0
    while True:

        for u in range(p):
            l[i][u] = b[i][u]-a[i][u]  

                                                                             
        if ( l[i][0] <= epsilon_metody ):  
            break
        if i+1 == L_zlotego_podzialu:                                              # tak w zapasie !!!!!!!!!!!!!!!!!!!
            break

        for u in range(p):
            v[i][u] = (a[i][u] + 0.5 * (3-math.sqrt(5)) * l[i][u])
            w[i][u] = (a[i][u] + 0.5 * (math.sqrt(5)-1) * l[i][u])

            zmienne4[u] = X[k][u] + v[i][u]*d[k][u]
            zmienne5[u] = X[k][u] + w[i][u]*d[k][u]

        if f( zmienne4 ) < f( zmienne5 ):  
            for u in range(p):
                a[i+1][u] = a[i][u]
                b[i+1][u] = w[i][u]
            i=i+1
        else:
            for u in range(p):
                a[i+1][u] = v[i][u]                
                b[i+1][u] = b[i][u]
            i=i+1
        #print("a = " + str(a[i]) + "oraz b = " + str(b[i]))
    for u in range(p):
        alfa[k][u] = ( (1/2) * ( a[i][u] + b[i][u]) )
        X[k+1][u]  = (X[k][u] + alfa[k][u]*d[k][u])
    print("alfa = " + str(alfa[k]))

    return X

def algorytm_bisekcji_glodensteina(poczatkowe_tau_r, beta, gradient_T, d, k, p, X, f):
    
    tau_l = 0
    tau_r = 0
    tau = 0
    tau_r = poczatkowe_tau_r
    czy_bylo = 'nie'
    
    nr_iteracji = 0

    zmienne1= [0] * p
    zmienne2= [0] * p
  
    #gradient_T = np.transpose(gradient[k])    
    #print('gradient_T = ' + str(gradient_T))
    p_kierunku = np.dot(gradient_T[k], d[k])
    print("p_kierunku = " +str(p_kierunku) + ' musi byc <0')
    while True:

        tau = 1/2 * (tau_l + tau_r)

        #print(nr_iteracji)
        #print('tau = ' + str(tau) )
        #print('taul = ' + str(tau_l) )
        #print('taur = ' + str(tau_r) )
        for i in range(p):
            zmienne1[i] = X[k][i]+ tau * d[k][i]
        #print(str(f(zmienne1))+' < '+str(f(X[k]) + (1-beta) * p_kierunku * tau))
       # print(str(f(zmienne1))+' > '+str(f(X[k]) + beta * p_kierunku * tau))
        #print('------------')
        if f(zmienne1) < f(X[k]) + (1-beta) * p_kierunku * tau:
            tau_l = tau

            czy_bylo = 'tak'
            #print('bylo 1')      
        elif f(zmienne1) > f(X[k]) + beta * p_kierunku * tau:
            tau_r = tau
            #print('bylo 2')
            czy_bylo = 'tak'
        if czy_bylo == 'nie' or nr_iteracji == 2000:   # POZNIEJ  USUNAC
            #print("zadzialalo")         
            break

        nr_iteracji = nr_iteracji + 1
        czy_bylo = 'nie'

    for i in range(p):
        X[k+1][i]  = (X[k][i] + tau*d[k][i])
    

    for i in range(p):
            zmienne2[i] = X[k][i]+ poczatkowe_tau_r * d[k][i]
    if f(zmienne2) < f(X[k]):
        print(str(X[k]) +' + '+ str(poczatkowe_tau_r) + ' * ' + str(d[k]))
        print(' Dobrze dobrane poczatkowe tau_R :  '+str(f(zmienne2))+' < '+str(f(X[k])))
    else:
        print(str(X[k]) +' + '+ str(poczatkowe_tau_r) + ' * ' + str(d[k]))
        print(' UWAGA !!! Poczatkowe tau_R nie spelnia warunku: f(x+tau_R*d) < f(x) : '+str(f(zmienne2))+' < '+str(f(X[k])))

    print('tau ostateczne =' + str(tau))
    print('d - kierunek =' + str(d[k]))
    print('X wczesniejsze =' + str(X[k]))
    print('X nowe =' + str(X[k+1]))
    return X
        


def rysowanie_wykresu(wykres_od, wykres_do, dokladnosc, F, X, wartosc_funkcji, k, p):                 
    
    if p==1:
        x = np.linspace(wykres_od,wykres_do,dokladnosc)
        F=str(F)
        Z = eval(F)

        plt.plot(x,Z,'r')

        for i in range(k+1):
            plt.plot(X[i],wartosc_funkcji[i], marker = '*')  
        plt.grid()
        plt.show()


    if p==2:         
        G = [] 
        H = [] 

        for i in range(k+1):
            if X[i][0]>wykres_od and X[i][0]<wykres_do and X[i][1]>wykres_od and X[i][1]<wykres_do:     # jezeli ten punkt znajduje sie w przedziale rysowania wykresu to przyjmij, jezeli nie to nie rysuj go, bo inaczej zepsuje wykres
                G.append(X[i][0]) 
                H.append(X[i][1])
       
        x = np.linspace(wykres_od,wykres_do,dokladnosc)              
        y = np.linspace(wykres_od,wykres_do,dokladnosc)          

        x,y = np.meshgrid(x, y)                        
        #F=str(F)
        print(F)
       
        # zeby rysowalo rownie z funkcje jak sin, cos, itd potrzeba uzyc biblioteki numpy, ponizszy kod dodaje ta biblioteke do odpowiednich funkcji
        j=0
        i=0
        Y=""
        while i != int(len(F)-1):
            if F[i].isalpha()==True and F[i+1].isalpha()==True: 
                j=j+1
            elif j!=0:
                F = F[0:i-j] + "np."+F[i-j:int(len(F))]    # dodaj Y zamiast F
                i=i+3
                j=0
            if i==int(len(F)-2) and j!=0:
                print('sprwadzenie1 - ' + str(F[0:i-j])+ ' - sprawdzenie2 - ' + str(F[i-j:int(len(F))]))
                F = F[0:i-j+1] + "np."+F[i-j+1:int(len(F))]   
                i=i+3
                j=0
            # if i==int(len(F)-1) nad F[i]:      
            i=i+1

                
        print(F)
        

        Z = eval(F)

        plt.figure(figsize=(10, 7))
        #plt.axes([0.1, 0.1, 0.8, 0.8])
        plt.plot(G,H,marker = '*') 
        ax1 = plt.contourf(x, y, Z, 30, alpha=.75, cmap=plt.cm.hot)     # wartosc po Z to liczbe warstwic na kazdy poziom, im wiecej tym wiecej wykreslonyc hwartwic i dokladniejszy wykres, ale uwazaj im wiecej tym dlzuej laduje wykres i wiecej wartswic nie oznacza ze wykres bedzie ladniejszy i czytelniejszy
        C = plt.contour(x, y, Z, 30, colors='black', linewidth=.5)
        plt.clabel(C, inline=1, fontsize=10)
        cbar = plt.colorbar(ax1, shrink=0.9)

        print("test" + str(G) + "oraz "  + str(H) )

        #plt.xticks([])    # wtedy znika nam oś x i y 
        #plt.yticks([])
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()



#def metoda_polaka_ribierry (punkt_poczatkowy, epsilon, epsilon1, epsilon2, epsilon_metody, m, L_zlotego_podzialu, L, p, F, znak_niewiadomej):
def metoda_polaka_ribierry (punkt_poczatkowy, epsilon, epsilon1, epsilon2, poczatkowe_tau_r, beta_gold, L, p, F, znak_niewiadomej):
    k=0
   
    # tworzenie tablic dwuwymiarowych
    X = [0] * L
    d = [0] * L
    beta = [0] * L
    gradient = [0] * L
    #tmp_grad = [0] * p
    hesjan = [0] * p
    wartosc_funkcji = [0] * L


    for i in range(L):
        X[i] = [0] * p
        d[i] = [0] * p
       # beta[i] = [0] * p                    # SPRAWDZIC CZY BETA JEST TYLKO WEKTOREM I USUSNAC !!!!!!!!!!!!!!
        gradient[i] = [0] * p

    for i in range(p):
        hesjan[i] = [0]* p


    for i in range(p):
        znak_niewiadomej[i] = sym.symbols(znak_niewiadomej[i])    
   
    # 1)

    for i in range(p):
        X[k][i] = float(punkt_poczatkowy[i])


    f = lambdify([znak_niewiadomej], F, 'numpy') 

    
    while True:
        # 2) 
        # wyliczamy wartosc funckji
     
        wartosc_funkcji[k] = f(X[k])   

        print('Krok: ' + str(k) + '; f(' + str(X[k]) + ') = ' + str(wartosc_funkcji[k]) )                                               #  USUN!!!!!!!!!!!!!!!!!!!
    
        # wyliczamy wzor na gradient   oraz    
        for i in range(p):
            gradient[k][i] = lambdify([znak_niewiadomej], (sym.diff(F,znak_niewiadomej[i])), 'numpy')     # (sym.diff(F,x)) = liczenie pochodnej funkcji F po niewiadomych x, nastepnie twozenie funckji gradient pod biblioteke numpy
        
         

        print(sym.diff(F,znak_niewiadomej[0])) 
        print(sym.diff(F,znak_niewiadomej[1])) 

        # wyliczamy wartosc gradientu
        for i in range(p):
            gradient[k][i] = gradient[k][i](X[k])   
        
        print('Wartosc gradientu=  ' + str(gradient[k]))
       
        # 3) 

   
        gradient_T = np.transpose(gradient[k]) # SPRAWDZIC BEZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        #print('Wartosc gradientu transponowanego=  ' + str(gradient_T))
   
        print("1) kryterium stopu = " + str(np.dot(gradient_T, gradient[k])))
        if np.dot(gradient_T, gradient[k]) <= epsilon:                  #??????? Sprawdzic tez czy na pewno dobrze to liczy
            print("STOP - kryterium stopu nr.1")
            kryterium_stopu = 'Epsylon - kryterium gradientowe ' 
            wartosc_kryterium_stopu = np.dot(gradient_T, gradient[k])
            break

        if k!=0: # gdyz na poczatku nie ma ani punktu ani wartosci do ktorej mozna by wyliczyc podane ponizej kryteria stopu

            print("2) kryterium stopu = " + str(math.sqrt( np.dot( np.transpose( np.subtract(X[k],X[k-1])) ,np.subtract(X[k],X[k-1]) ) ) ))       #math.sqrt( np.dot( np.transpose(X[k]) ,X[k] ) ) - math.sqrt(np.dot( np.transpose(X[k-1]) ,X[k-1] ))
            if math.sqrt( np.dot( np.transpose( np.subtract(X[k],X[k-1])) ,np.subtract(X[k],X[k-1]) ) )   <= epsilon1:        #??????? WSZYSTKO POD JEDEN PIERWIASTEK I POZNIZEJ TAKZE I SPRAWDZ WSZYSTKIE WARUNKI CZY DOBRZE ZROBIONE
                print("STOP - kryterium stopu nr.2")
                kryterium_stopu = 'Epsylon1 - kryterium  normy argumentów' 
                wartosc_kryterium_stopu = math.sqrt( np.dot( np.transpose(np.subtract(X[k],X[k-1])) ,np.subtract(X[k],X[k-1])) ) 
                break
            print("3) kryterium stopu = " + str (abs(wartosc_funkcji[k]-wartosc_funkcji[k-1]))) #(math.sqrt( np.dot( np.transpose( np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1])) ,np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1])) )))                                                                                                  #math.sqrt( np.dot( np.transpose(wartosc_funkcji[k]) ,wartosc_funkcji[k] ) ) - math.sqrt(np.dot( np.transpose(wartosc_funkcji[k-1]) ,wartosc_funkcji[k-1] ))
            if abs(wartosc_funkcji[k]-wartosc_funkcji[k-1])<= epsilon2:  #math.sqrt( np.dot( np.transpose( np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1])) ,np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1])) )  <= epsilon2:                 #???????
                print("STOP - kryterium stopu nr.3")
                kryterium_stopu = 'Epsylon2 - kryterium  wartosci funkcji' 
                wartosc_kryterium_stopu = abs(wartosc_funkcji[k]-wartosc_funkcji[k-1]) #math.sqrt( np.dot( np.transpose(np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1])) ,np.subtract(wartosc_funkcji[k],wartosc_funkcji[k-1]) ) )
                break
       
        if k+1 == L:                                              
            print("STOP - kryterium maksimum iteracji")
            kryterium_stopu = 'L - kryterium maksimum iteracji ???' 
            wartosc_kryterium_stopu = k
            break
        
       

        #4) Kierunek poszukiwan
            
        d = kierunek_poszukiwan(k, gradient,p,beta,d)

        #5)  algorytm zlotego podzialu
        #X = metoda_zlotego_podzialu(epsilon_metody, m, L_zlotego_podzialu, X, d, k,p,f)

        #5)  algorytm goldensteina
        
        X = algorytm_bisekcji_glodensteina(poczatkowe_tau_r, beta_gold, gradient, d, k, p, X, f)
        #for i in range(p): ##########################################
        #X[2][i]  = (X[1][i] + 0.008*d[1][i]) 
        k=k+1


    
    # Wyliczamy hesjan by sprawdzic czy punkt na pewno znajduje sie w min funkcji (gdyz moze zarowno znajdowac sie na wierzcholku lub na przegieciu funkcji)

    for i in range(p):
        tmp_grad = (sym.diff(F,znak_niewiadomej[i]))
        for l in range(p):
            tmp_hesjan = sym.diff(tmp_grad,znak_niewiadomej[l])
            print(tmp_hesjan)
            hesjan[i][l] = lambdify([znak_niewiadomej], tmp_hesjan  , 'numpy') 
    for i in range(p):
        for l in range(p):
            hesjan[i][l] = hesjan[i][l](X[k]) 
            print(hesjan[i][l])

    if np.linalg.det(hesjan) > 0:
        print("Tutaj HESJAN, jest ok :)")
    else:
        print("UWAGA!!!  Hesjan NIE jest dodatnio określony, znaleziony punkt nie jest minimum")
    print('Wartosc wyznacznika hesjanu = ' + str(np.linalg.det(hesjan)))
    print('punkt koncowy X = ' + str(X[k]))
    print('Wartosc koncowa funkcji = ' + str(wartosc_funkcji[k]))                                                                                #USUN!!!!!!!!!!!!!!!!!!!


    return X, wartosc_funkcji, k,  kryterium_stopu, wartosc_kryterium_stopu



'''
1) Wybierz punk startowy x0, podaj funkcje f(x), przypisz punk startowy x[k]=x[0]
2) Oblicz wartosc funkcji f(x[k]) oraz jezeli jest to wymagane to jej gradient rowniez 
3) Zbadaj przyjęte kryteria
      Jeżeli kryterium spełnione zakończ wykonywanie algorytmu - uzyskano rozwiazanie optymalne x[k] i potymalna wartosc funkcji celu f(x[k])
      Jeżeli nie to przejdz do 4)
4) Wyznacz ustalony kierunek poszukiwań: d[k] (metoda polaka-ribiery)
5) Wykonaj minimalizację kierunkową wybraną metodą (metoda minimum kierunku)
6) Podstaw x[k]=x[k+1]   oraz k=k+1
7) Przejdz do 2)


minimalizacja kierunku / metoda minimum kierunku - metoda złotego podziału

http://smurf.mimuw.edu.pl/external_slides/MO/MO_6_metody_rozwiazywania/MO_6_metody_rozwiazywania.html

Inicjalizacja:

Podaj pcozatkowa dlugosc kroku m>0
Wybierz dokladnosc lokalizacji minimum epsilon>0

1) Oblicz f(x[k] + m*d[k])
2) Jezeli f(x[k] + m*d[k]) >= f(x[k]) podstaw a[0] = 0, b[0] = m i idź do 5)
   Jezeli nie, podstaw j=0 alfa[j]=0 i idz do 3)
3) Podstaw alfa[j+1]=alfa[j] + m i oblicz f(x[k]+ alfa[j+1] * d[k])
4) Jeżeli f(x[k]+ alfa[j+1] * d[k]) >= f(x[k]+ alfa[j] * d[k]) podstaw a[0]=alfa[j-1], b[0]=alfa[j+1] i idz do 5)
   Jezeli nie podstaw j=j+1 i idz do 3)
5) Podstaw i=0
6) Podstaw l[i]=b-a
7) Jezeli L[i]<espilon idź do 10).
   Jezeli nie, idz do 8
8) Podstaw v[i]=a[i] + 1/2 * (3-sqrt(5)) * l[i] oraz w[i]=a[i] + 1/2 * (sqrt(5)-1) * l[i]
   Oblicz f( x[k] + v[i]*d[k] ) oraz f( x[k] + w[i]*d[k] )
9) Jezeli f( x[k] + v[i]*d[k] ) < f( x[k] + w[i]*d[k] , podstaw a[i+1]=a[i], b[i+1]=w[i], i=i+1 i idz do 6)
   Jezeli nie to podstaw a[i+1]=v[i], b[i+1]=b[i] , i=i+1 i idz do 6)
10) Przyjmij za rozwiazanie alfa[k]= 1/2  * ( a[i]*b[i] ) i stop.



Jak to dziala na logike:
w punktach 1)-4) szukamy przedzialu w ktorym moze znajdowac sie minimum lokalne. W tym celu sprawdzamy
czy punkt kolejny ma mniejsza wartosc od poprzedniego (kolejny punkt wyznaczamy na podstawie poczatkowej 
dlugosci "m" oraz kierunku poszukiwan "d"), jezeli TAK to znalezlismy przedzial z min. lokal., jak NIE (kolejny 
punkt jest mniejszy od poprzedniego, czyli funkcja maleje) to przesuwamy te dwa punkty o pocz. dlug. m
i ponownie sprawdzamy. I tak w kolko az nie znajdziemy minimum.
Kiedy znajdzeimy minimum przechodzimy do zawezania go ( punkty 6)-9) ), az nie osiagniemy odpowiednio
dokladna wartosc.


Do testow iteracji wez x^3 - x + 1
'''






