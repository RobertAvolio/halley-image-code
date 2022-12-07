import random

def main():

    random.seed(30)
    
    loops = input("How many files to generate?\n")
    
    for i in range(int(loops)):
        
        name = "testfiles/test" + str(i) + ".txt"
        f = open(name, "w")
        
        helio = random.randrange(0, 50, 1)
        f.write(helio "\n")
        didist = random.randrange(0, 50, 1)
        f.write(didist + "\n")
        scale = random.randrange(, , )
        f.write(scale + "\n")
        cenlon = random.randrange(-90, 90, 1)
        f.write(cenlon + "\n")
        cenlat = random.randrange(-90, 90, 1)
        f.write(cenlat + "\n")
        raddeg = random.randrange(0, 90, .1)
        f.write(raddeg + "\n")
        sigjet = random.randrange(20, 50, 3)
        f.write(sigjet + "\n")
        angra = random.randrange(-90, 90, 1)
        f.write(angra + "\n")
        angdec = random.randrange(-90, 90, 1)
        f.write(angdec + "\n")
        solra = random.randrange(-90, 90, 1)
        f.write(solra + "\n")
        soldec = random.randrange(-90, 90, 1)
        f.write(soldec + "\n")
        earra = random.randrange(-90, 90, 1)
        f.write(earra + "\n")
        eardec = random.randrange(-90, 90, 1)
        f.write(eardec + "\n")
        numap = random.randrange(0, 1, 1)
        f.write(numap + "\n")
        for i in range(numap):
            rad = random.randrange(0, 1, 1)
            f.write(rad + "\n")
            
        angthe = random.randrange(0, 90, 1)
        f.write(angthe + "\n")
        
        ratpsi = random.randrange(0, 1, 1)
        f.write(ratpsi + "\n")
        ratphi = random.randrange(0, 0, 1)
        f.write(ratphi + "\n")
        date = random.randrange(0, 0, 1)
        f.write(date + "\n")
        gpsi = random.randrange(0, 0, 1)
        f.write(gpsi + "\n")
        gphi = random.randrange(0, 0, 1)
        f.write(gphi + "\n")
        radg = random.randrange(0, 0, 1)
        f.write(radg + "\n")
        radr = random.randrange(0, 0, 1)
        f.write(radr + "\n")
        velgr = random.randrange(0, 0.5, .1)
        f.write(velgr + "\n")
        velr = random.randrange(, , 1)
        f.write(velr + "\n")
        lifgra = random.randrange(, , 1)
        f.write(lifegra + "\n")
        lifrad = random.randrange(, , 1)
        f.write(lifrad + "\n")
        timbin = random.randrange(, , 1)
        f.write(timbin + "\n")
        totime = random.randrange(, , 1)
        f.write(totime + "\n")
        nregion = random.randrange(, , 1)
        f.write(nregion + "\n")
        binnum = random.randrange(, , 1)
        f.write(binnum + "\n")
        bg = random.randrange(, , 1)
        f.write(bg + "\n")
        f.close

if __name__ == "__main__":
    main()
