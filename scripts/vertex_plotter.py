import sys
sys.path.insert(1, '..')

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np
import ROOT
from functions.other_functions import io_parse_arguments


def fill_gveto_om(gveto_om_list: list, side_num: int, end_num: int, y_vertex_pos: float):
    s = side_num
    t = end_num
    z = y_vertex_pos
    if z > 0:
        if 0 < z < 330:
            om = 8
        elif 330 < z < 630:
            om = 9
        elif 630 < z < 930:
            om = 10
        elif 930 < z < 1230:
            om = 11
        elif 1230 < z < 1530:
            om = 12
        elif 1530 < z < 1830:
            om = 13
        elif 1830 < z < 2130:
            om = 14
        elif z > 2130:
            om = 15
    elif z < 0:
        r_z = z
        z = -z
        if 0 < z < 330:
            om = 7
        elif 330 < z < 630:
            om = 6
        elif 630 < z < 930:
            om = 5
        elif 930 < z < 1230:
            om = 4
        elif 1230 < z < 1530:
            om = 3
        elif 1530 < z < 1830:
            om = 2
        elif 1830 < z < 2130:
            om = 1
        elif z > 2130:
            om = 0
    gveto_om_list[s][t][om] += 1


def get_xwall_col_num(x_vertex_pos: int):
    x = abs(x_vertex_pos)
    if x < 230:
        col_num = 0
    else:
        col_num =1
    return col_num


def fill_xwall_om(xwall_om_list: list, side_num: int, end_num: int, column_num: int, z_vertex_pos: float):
    z = z_vertex_pos
    s = side_num
    e = end_num
    c = column_num

    if z > 0:
        if z < 210:
            om_num = 8
        elif 210 < z < 420:
            om_num = 9
        elif 420 < z < 630:
            om_num = 10
        elif 630 < z < 840:
            om_num = 11
        elif 840 < z < 1050:
            om_num = 12
        elif 1050 < z < 1260:
            om_num = 13
        elif 1260 < z < 1470:
            om_num = 14
        elif z > 1470:
            om_num = 15
    else:
        r_z = z
        z = -z
        if z < 210:
            om_num = 7
        elif 210 < z < 420:
            om_num = 6
        elif 420 < z < 630:
            om_num = 5
        elif 630 < z < 840:
            om_num = 4
        elif 840 < z < 1050:
            om_num = 3
        elif 1050 < z < 1260:
            om_num = 2
        elif 1260 < z < 1470:
            om_num = 1
        elif z > 1470:
            om_num = 0
    xwall_om_list[s][e][c][om_num] += 1


def fill_mwall_om(mwall_om_list: list, side_num: int, y_vertex_pos: float, z_vertex_pos: float):
    y = y_vertex_pos
    z = z_vertex_pos

    # Isolate column
    if y < -2350:
        col = 0
    elif -2350 < y < -2080:
        col = 1
    elif -2080 < y < -1820:
        col = 2
    elif -1820 < y < -1560:
        col = 3
    elif -1560 < y < -1300:
        col = 4
    elif -1300 < y < -1040:
        col = 5
    elif -1040 < y < -780:
        col = 6
    elif -780 < y < -530:
        col = 7
    elif -530 < y < -260:
        col = 8
    elif -260 < y < 0:
        col = 9
    elif y > 2350:
        col = 19
    elif 2350 > y > 2080:
        col = 18
    elif 2080 > y > 1820:
        col = 17
    elif 1820 > y > 1560:
        col = 16
    elif 1560 > y > 1300:
        col = 15
    elif 1300 > y > 1040:
        col = 14
    elif 1040 > y > 780:
        col = 13
    elif 780 > y > 530:
        col = 12
    elif 530 > y > 260:
        col = 11
    elif 260 > y > 0:
        col = 10
    if z < -1440:
        row = 0
    elif -1440 < z < -1150:
        row = 1
    elif -1150 < z < -890:
        row = 2
    elif -890 < z < -640:
        row = 3
    elif -640 < z < -380:
        row = 4
    elif -380 < z < -120:
        row = 5
    elif -120 < z < 120:
        row = 6
    elif 120 < z < 380:
        row = 7
    elif 380 < z < 640:
        row = 8
    elif 640 < z < 890:
        row = 9
    elif 890 < z < 1150:
        row = 10
    elif 1150 < z < 1440:
        row = 11
    elif 1440 < z:
        row = 12

    if row == 0 or row == 12:
        n = 5
    else:
        n = 8

    mwall_om_list[side_num][col][row] += 1

    return n


def main():
    args = io_parse_arguments()
    input_data_filename = args.i
    output_filename = args.o

    name = input_data_filename.split("/")[-1]
    name = name.split(".root")[0]

    gveto_num = 0
    xwall_num = 0
    mw_8i_num = 0
    mw_5i_num = 0
    mwall_num = 0

    x_array_gveto = []
    y_array_gveto = []
    z_array_gveto = []

    x_array_xwall = []
    y_array_xwall = []
    z_array_xwall = []

    x_array_mw_8i = []
    y_array_mw_8i = []
    z_array_mw_8i = []

    x_array_mw_5i = []
    y_array_mw_5i = []
    z_array_mw_5i = []

    gveto_om_list = [[[0 for i in range(16)], [0 for i in range(16)]], [[0 for i in range(16)],  [0 for i in range(16)]]]
    xwall_om_list = [[[[0 for i in range(16)], [0 for i in range(16)]], [[0 for i in range(16)], [0 for i in range(16)]]], [[[0 for i in range(16)], [0 for i in range(16)]], [[0 for i in range(16)], [0 for i in range(16)]]]]
    mwall_om_list = [[[0 for i in range(13)] for j in range(20)], [[0 for i in range(13)] for j in range(20)]]

    file = ROOT.TFile(input_data_filename, "READ")
    output = ROOT.TFile(output_filename,"RECREATE")

    xy_hist = ROOT.TH2I("xy", "xy", 8000, -4000, 4000, 8000, -4000, 4000)
    yz_hist = ROOT.TH2I("yz", "yz", 8000, -4000, 4000, 8000, -4000, 4000)
    xz_hist = ROOT.TH2I("xz", "xz", 8000, -4000, 4000, 8000, -4000, 4000)

    gveto_event_rate = ROOT.TH1I("gveto_rate", "gveto_rate", 1000000, 0, 1000000)
    xwall_event_rate = ROOT.TH1I("xwall_rate", "xwall_rate", 1000000, 0, 1000000)
    mwall_event_rate_8inch = ROOT.TH1I("mwall_rate_8inch", "mwall_rate_8inch", 1000000, 0, 1000000)
    mwall_event_rate_5inch = ROOT.TH1I("mwall_rate_5inch", "mwall_rate_5inch", 1000000, 0, 1000000)

    om_2d_hist = ROOT.TH2I(name, name, 24, 0, 24, 36, 0, 36)

    for event in file.T:
        x = event.x
        y = event.y
        z = event.z

        xy_hist.Fill(x, y)
        yz_hist.Fill(y, z)
        xz_hist.Fill(x, z)

        # Cuts on xy

        # Isolate XW and MW OMs
        if (-1700 < z < 1700):
            # XW
            if (-550 < x < 550):
                x_array_xwall.append(x)
                y_array_xwall.append(y)
                z_array_xwall.append(z)
                xwall_num += 1
                # End 1
                if y > 0:
                    e = 1
                    # Side 1
                    if x < 0:
                        s = 1
                    # Side 0
                    else:
                        s = 0

                # End  0
                else:
                    e = 0
                    # Side 1
                    if x < 0:
                        s = 1
                    # Side 0
                    else:
                        s = 0
                c = get_xwall_col_num(x)
                fill_xwall_om(xwall_om_list, s, e, c, z)
            # MW
            else:
                mwall_num += 1
                # Side 1
                if x > 0:
                    s = 1
                # Side 0
                else:
                    s = 0
                type = fill_mwall_om(mwall_om_list, s, y, z)
                if type == 5:
                    mw_5i_num += 1
                    x_array_mw_5i.append(x)
                    y_array_mw_5i.append(y)
                    z_array_mw_5i.append(z)
                else:
                    mw_8i_num += 1
                    x_array_mw_8i.append(x)
                    y_array_mw_8i.append(y)
                    z_array_mw_8i.append(z)
        # The rest are the GV OMs
        else:
            gveto_num += 1
            x_array_gveto.append(x)
            y_array_gveto.append(y)
            z_array_gveto.append(z)
            # Top
            if z > 0:
                t = 1
                # Side 1
                if x > 0:
                    s = 1
                # Side 0
                else:
                    s = 0
            # Bottom
            else:
                t = 0
                # Side 1
                if x > 0:
                    s = 1
                # Side 0
                else:
                    s = 0
            fill_gveto_om(gveto_om_list, s, t, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(np.array(x_array_mw_5i), np.array(y_array_mw_5i), np.array(z_array_mw_5i), c='g', marker='.', label='MW 5\"')
    ax.scatter(np.array(x_array_mw_8i), np.array(y_array_mw_8i), np.array(z_array_mw_8i), c='k', marker='.', label='MW 8\"')
    ax.scatter(np.array(x_array_gveto), np.array(y_array_gveto), np.array(z_array_gveto), c='r', marker='.', label='GVeto')
    ax.scatter(np.array(x_array_xwall), np.array(y_array_xwall), np.array(z_array_xwall), c='b', marker='.', label='XWall')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.legend()
    plt.savefig('/Users/willquinn/Desktop/pmt_glass_sim_test/' + name)

    plt.close(fig)

    gveto_oms = []
    xwall_oms = []
    mw_5i_oms = []
    mw_8i_oms = []

    for i_side in range(2):
        if i_side == 0:
            s = -1
        else:
            s = 1
        for i_end in range(2):
            for i_col in range(2):
                for i_om in range(len(xwall_om_list[i_side][i_end][i_col])):
                    if i_end == 1:
                        x = 23 - i_col
                        if i_side == 0:
                            y = i_om + 1
                        else:
                            y = i_om + 19
                    else:
                        x = 0 + i_col
                        if i_side == 0:
                            y = i_om + 1
                        else:
                            y = i_om + 19
                    xwall_event_rate.Fill(xwall_om_list[i_side][i_end][i_col][i_om])
                    xwall_oms.append(xwall_om_list[i_side][i_end][i_col][i_om])
                    for i in range(xwall_om_list[i_side][i_end][i_col][i_om]):
                        om_2d_hist.Fill(x, y)

            for i_om in range(len(gveto_om_list[i_side][i_end])):
                if i_side == 0:
                    if i_end == 1:
                        y = 17
                    else:
                        y = 0
                else:
                    if i_end == 1:
                        y = 35
                    else:
                        y = 18
                x = i_om + 4
                gveto_event_rate.Fill(gveto_om_list[i_side][i_end][i_om])
                gveto_oms.append(gveto_om_list[i_side][i_end][i_om])
                for i in range(gveto_om_list[i_side][i_end][i_om]):
                    om_2d_hist.Fill(x, y)

        for i_col in range(len(mwall_om_list[i_side])):
            for i_row in range(len(mwall_om_list[i_side][i_col])):
                if i_side == 0:
                    y = 2 + i_row
                else:
                    y = 20 + i_row
                x = 2 + i_col
                if i_row == 12 or i_row == 0:
                    mwall_event_rate_5inch.Fill(mwall_om_list[i_side][i_col][i_row])
                    mw_5i_oms.append(mwall_om_list[i_side][i_col][i_row])
                else:
                    mwall_event_rate_8inch.Fill(mwall_om_list[i_side][i_col][i_row])
                    mw_8i_oms.append(mwall_om_list[i_side][i_col][i_row])
                for i in range(mwall_om_list[i_side][i_col][i_row]):
                    om_2d_hist.Fill(x, y)

    gveto_oms_array = np.array(gveto_oms)
    g_x = []
    for i in range(gveto_oms_array.size):
        g_x.append(i)
    g_x = np.array(g_x)

    xwall_oms_array = np.array(xwall_oms)
    x_x = []
    for i in range(xwall_oms_array.size):
        x_x.append(i)
    x_x = np.array(x_x)

    mw_5i_oms_array = np.array(mw_5i_oms)
    m_x_5 = []
    for i in range(mw_5i_oms_array.size):
        m_x_5.append(i)
    m_x_5 = np.array(m_x_5)

    mw_8i_oms_array = np.array(mw_8i_oms)
    m_x_8 = []
    for i in range(mw_8i_oms_array.size):
        m_x_8.append(i)
    m_x_8 = np.array(m_x_8)

    g_test = 0
    x_test = 0
    m_test_5 = 0
    m_test_8 = 0

    if gveto_num > 0:
        g_test = 1
    if xwall_num > 0:
        x_test = 1
    if mw_5i_num > 0:
        m_test_5 = 1
    if mw_8i_num > 0:
        m_test_8 = 1

    num_gveto_expected = (gveto_num + xwall_num + mw_8i_num + mw_5i_num)/(g_test*64 + x_test*128 + m_test_5*80 + m_test_8*(189.444/101.832)*440)
    g_sig = np.sqrt(num_gveto_expected)
    num_xwall_expected = (gveto_num + xwall_num + mw_8i_num + mw_5i_num)/(g_test*64 + x_test*128 + m_test_5*80 + m_test_8*(189.444/101.832)*440)
    x_sig = np.sqrt(num_xwall_expected)
    num_mw_5i_expected = (gveto_num + xwall_num + mw_8i_num + mw_5i_num)/(g_test*64 + x_test*128 + m_test_5*80 + m_test_8*(189.444/101.832)*440)
    m_sig_5 = np.sqrt(num_mw_5i_expected)
    num_mw_8i_expected = (gveto_num + xwall_num + mw_8i_num + mw_5i_num)/(g_test* (101.832/189.444) *64 + x_test*(101.832/189.444)*128 + m_test_5*(101.832/189.444)*80 + m_test_8*440)
    m_sig_8 = np.sqrt(num_mw_8i_expected)

    if g_test == 1:
        plt.errorbar(g_x, gveto_oms_array, yerr=np.sqrt(gveto_oms_array), fmt='k.', ecolor='k')
        plt.plot([], [], 'k.', label="GVeto Num of Vertexes")
        plt.plot(g_x, np.ones_like(g_x)*num_gveto_expected, "r-", label="Expected = {:.3f}".format(num_gveto_expected))
        plt.plot(g_x, np.ones_like(g_x) * (num_gveto_expected + g_sig), "b-", label="+/- 1 $\sigma$")
        plt.plot(g_x, np.ones_like(g_x) * (num_gveto_expected - g_sig), "b-")
        plt.ylabel("Number of vertexes")
        plt.legend()
        plt.savefig('/Users/willquinn/Desktop/pmt_glass_sim_test/' + name + "_gvet")
        plt.close()

    if x_test == 1:
        plt.errorbar(x_x, xwall_oms_array, yerr=np.sqrt(xwall_oms_array), fmt='k.', ecolor='k')
        plt.plot([], [], 'k.', label="XWall Num of Vertexes")
        plt.plot(x_x, np.ones_like(x_x) * num_xwall_expected, "r-", label="Expected = {:.3f}".format(num_xwall_expected))
        plt.plot(x_x, np.ones_like(x_x) * (num_xwall_expected + x_sig), "b-", label="+/- 1 $\sigma$")
        plt.plot(x_x, np.ones_like(x_x) * (num_xwall_expected - x_sig), "b-")
        plt.ylabel("Number of vertexes")
        plt.legend()
        plt.savefig('/Users/willquinn/Desktop/pmt_glass_sim_test/' + name + "_xwal")
        plt.close()

    if m_test_5 == 1:
        plt.errorbar(m_x_5, mw_5i_oms_array, yerr=np.sqrt(mw_5i_oms_array), fmt='k.', ecolor='k')
        plt.plot([], [], 'k.', label="MWall 5\" Num of Vertexes")
        plt.plot(m_x_5, np.ones_like(m_x_5) * num_mw_5i_expected, "r-", label="Expected = {:.3f}".format(num_mw_5i_expected))
        plt.plot(m_x_5, np.ones_like(m_x_5) * (num_mw_5i_expected + m_sig_5), "b-", label="+/- 1 $\sigma$")
        plt.plot(m_x_5, np.ones_like(m_x_5) * (num_mw_5i_expected - m_sig_5), "b-")
        plt.ylabel("Number of vertexes")
        plt.legend()
        plt.savefig('/Users/willquinn/Desktop/pmt_glass_sim_test/' + name + "_mw5i")
        plt.close()

    if m_test_8 == 1:
        plt.errorbar(m_x_8, mw_8i_oms_array, yerr=np.sqrt(mw_8i_oms_array), fmt='k.', ecolor='k')
        plt.plot([], [], 'k.', label="MWall 8\" Num of Vertexes")
        plt.plot(m_x_8, np.ones_like(m_x_8) * num_mw_8i_expected, "r-", label="Expected = {:.3f}".format(num_mw_8i_expected))
        plt.plot(m_x_8, np.ones_like(m_x_8) * (num_mw_8i_expected + m_sig_8), "b-", label="+/- 1 $\sigma$")
        plt.plot(m_x_8, np.ones_like(m_x_8) * (num_mw_8i_expected - m_sig_8), "b-")
        plt.ylabel("Number of vertexes")
        plt.legend()
        plt.savefig('/Users/willquinn/Desktop/pmt_glass_sim_test/' + name + "_mw8i")
        plt.close()

    c1 = ROOT.TCanvas(name)
    c1.cd()
    c1.SetGrid()
    om_2d_hist.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)
    c1.SaveAs("~/Desktop/pmt_glass_sim_test/" + name + ".pdf")

    print(">>> Number of simulated vertexes:          ", gveto_num + xwall_num + mw_8i_num + mw_5i_num)
    print(">>> Percentage of GVeto vertexes:          ", gveto_num/(gveto_num + xwall_num + mw_8i_num + mw_5i_num) * 100)
    print(">>> Mean number of vertexes per GVeto OM:  ", gveto_event_rate.GetMean())
    print(">>> StdDev of vertexes over GVeto OMs:     ", gveto_event_rate.GetStdDev())
    if gveto_event_rate.GetMean() == 0:
        pass
    else:
        print(">>> Resolution of GVeto vertexes:          ", gveto_event_rate.GetStdDev() / gveto_event_rate.GetMean() * 100)
    print(">>> Percentage of XWall vertexes:          ", xwall_num/(gveto_num + xwall_num + mw_8i_num + mw_5i_num) * 100)
    print(">>> Mean number of vertexes per XWall OM:  ", xwall_event_rate.GetMean())
    print(">>> StdDev of vertexes over XWall OMs:     ", xwall_event_rate.GetStdDev())
    if xwall_event_rate.GetMean() == 0:
        pass
    else:
        print(">>> Resolution of XWall vertexes:          ", xwall_event_rate.GetStdDev() / xwall_event_rate.GetMean() * 100)
    print(">>> Percentage of MWall 8\" vertexes:       ", mw_8i_num / (gveto_num + xwall_num + mw_8i_num + mw_5i_num) * 100)
    print(">>> Mean number of vertexes per MW 8\" OM:  ", mwall_event_rate_8inch.GetMean())
    print(">>> StdDev of vertexes over MW 8\" OMs:     ", mwall_event_rate_8inch.GetStdDev())
    if mwall_event_rate_8inch.GetMean() == 0:
        pass
    else:
        print(">>> Resolution of MW \8"" vertexes:          ",
              mwall_event_rate_8inch.GetStdDev() / mwall_event_rate_8inch.GetMean() * 100)
    print(">>> Percentage of MWall 5\" vertexes:       ", mw_5i_num / (gveto_num + xwall_num + mw_8i_num + mw_5i_num) * 100)
    print(">>> Mean number of vertexes per MW 5\" OM:  ", mwall_event_rate_5inch.GetMean())
    print(">>> StdDev of vertexes over MW 5\" OMs:     ", mwall_event_rate_5inch.GetStdDev())
    if mwall_event_rate_5inch.GetMean() == 0:
        pass
    else:
        print(">>> Resolution of MW 5\" vertexes:          ",
              mwall_event_rate_8inch.GetStdDev() / mwall_event_rate_5inch.GetMean() * 100)

    file.Close()
    output.cd()
    output.Write()
    output.Close()


if __name__ == '__main__':
    main()
