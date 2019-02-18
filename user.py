import hcsc
import photon_in
import test
import nu


def main():
    scell = hcsc.hcsc()
    ph_in = photon_in.Photon_in()
    ph_in.c = 10
    scell.shine(ph_in)
    scell.display_attributes()
    print("Jout:{:.3f}(A/m^2)".format(scell.Jouthc(0.2/2*nu.eV)[0]))
    print("Pout:{:.3f}(W/m^2)".format(scell.Pouthc(0.2/2*nu.eV)[0]))
    # print("MaxPout:{:.3f}(W/m^2)".format(scell.maxPouthc()[0]))

    test.IV(scell)
    test.PV(scell)
    test.maxPouthc_opt_ESC(scell, 50)
    return


if __name__ == '__main__':
    main()
