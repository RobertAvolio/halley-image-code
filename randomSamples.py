import random

def main():

    random.seed(30)

    loops = input("How many files to generate?\n")
    
    for i in range(int(loops)):
        
        name = "testfiles/test" + str(i) + ".txt"
        f = open(name, "w")
        
        helio = random.uniform(0, 50)
        f.write(str(helio) + "\n")
        didist = random.uniform(0, 50)
        f.write(str(didist) + "\n")
        scale = 0.3
        f.write(str(scale) + "\n")
        cenlon = random.uniform(-90, 90)
        f.write(str(cenlon) + "\n")
        cenlat = random.uniform(-90, 90)
        f.write(str(cenlat) + "\n")
        raddeg = random.uniform(0.1, 15)
        f.write(str(raddeg) + "\n")
        sigjet = random.uniform(0.1, 30)
        f.write(str(sigjet) + "\n")
        angra = random.uniform(-90, 90)
        f.write(str(angra) + "\n")
        angdec = random.uniform(-90, 90)
        f.write(str(angdec) + "\n")
        solra = random.uniform(-90, 90)
        f.write(str(solra) + "\n")
        soldec = random.uniform(-90, 90)
        f.write(str(soldec) + "\n")
        earra = random.uniform(-90, 90)
        f.write(str(earra) + "\n")
        eardec = random.uniform(-90, 90)
        f.write(str(eardec) + "\n")
        numap = random.randrange(0, 1, 1)
        # uses randrange so numap is an int
        f.write(str(numap) + "\n")
        for i in range(numap):
            rad = random.uniform(0, 1)
            f.write(str(rad) + "\n")
            
        angthe = random.uniform(0, 90)
        f.write(str(angthe) + "\n")
        
        ratpsi = random.uniform(0, 3024000)
        f.write(str(ratpsi) + "\n")
        ratphi = 0.5
        f.write(str(ratphi) + "\n")
        date = random.uniform(86400, 864000)
        f.write(str(date) + "\n")
        gpsi = 0
        f.write(str(gpsi) + "\n")
        gphi = 0
        f.write(str(gphi) + "\n")
        radg = 0
        f.write(str(radg) + "\n")
        radr = 0
        f.write(str(radr) + "\n")
        velgr = 0.5
        f.write(str(velgr) + "\n")
        velr = 0.00001
        f.write(str(velr) + "\n")
        lifgra = 0.00001
        f.write(str(lifgra) + "\n")
        lifrad = 100000
        f.write(str(lifrad) + "\n")
        timbin = random.uniform(1800, 86400)
        f.write(str(timbin) + "\n")
        totime = random.uniform(86400, 864000)
        while date >= totime:
            date = random.uniform(86400, 864000-1)
        f.write(str(totime) + "\n")
        # uses randrange so nregion/binnum are ints, 10001 not included so 10000 is the max
        nregion = random.randrange(100, 10001, 1)
        f.write(str(nregion) + "\n")
        binnum = random.randrange(100, 10000, 1)
        f.write(str(binnum) + "\n")
        bg = random.uniform(0, 16384)
        f.write(str(bg) + "\n")
        f.close()

if __name__ == "__main__":
    main()
