import sys

sys.path.insert(1, '../')

import numpy as np
import ROOT
from calomissioning.sndisplay import *
import matplotlib.pyplot as plt
from functions.other_functions import om_id_string, io_parse_arguments, inner_product


def sweep_example(om_text: str, template: np.array, waveform: np.array):
    sweep_start = 400 * 0.64
    mf_output = []

    x = np.array([i for i in range(waveform.size)]) / 0.64
    mf_x = []

    i = 0
    j = 0

    while i < waveform.size - template.size:
        mf_x.append(i / 0.64)
        test_x = np.array([j / 0.64 for j in range(i, i + template.size)])
        test = waveform[i:i + template.size]
        test_norm = np.sqrt(np.dot(test, test))
        shape = np.dot(template, test) / test_norm
        amplitude = np.dot(template, test)
        mf_output.append(shape)

        fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
        gs = fig.add_gridspec(3, hspace=0.3)
        axs = gs.subplots()
        axs[0].plot(x, waveform, "k")
        axs[0].set_title(om_text)
        axs[0].set_ylim([-900, 100])
        axs[0].set_xlim([x[0], x[-1]])
        axs[0].fill([test_x[0], test_x[0], test_x[-1], test_x[-1]], [100, -100, -100, 100], color='g', alpha=0.2)
        axs[0].set_ylabel("adc /mV")
        ax_0 = axs[0].twinx()
        ax_0.axvline(test_x[0], -4000, 4000, ls='--', color='g', label='sweep')
        ax_0.plot(test_x, template * test_norm, 'r', label='template')
        ax_0.plot(0, "k", label='waveform')
        ax_0.set_ylim([-900, 100])
        ax_0.set_ylabel('a.u')
        plt.legend(loc='lower right')

        axs[1].plot(test_x, test, "k")
        axs[1].set_ylabel("adc /mV")
        axs[1].set_xlim(test_x[0], test_x[-1])
        axs[1].set_ylim(np.min(test) + np.min(test) * 0.1 - 10,
                        -np.min(test) * 0.1 + 10)
        axs[1].grid()
        axs[1].set_title(f'mf shape = {shape}')
        ax_1 = axs[1].twinx()
        ax_1.set_ylim(np.min(test) + np.min(test) * 0.1 - 10,
                      -np.min(test) * 0.1 + 10)
        ax_1.set_ylabel('a.u')
        ax_1.plot(test_x, template * test_norm, "r")

        axs[2].plot(mf_x, mf_output, 'g.')
        axs[2].set_xlabel("timestamp /ns")
        axs[2].set_ylabel("shape")
        axs[0].grid()
        axs[2].grid()
        axs[2].set_title("convolution")
        plt.savefig(f"plots/sweep_plots/{j}.png")
        '''plt.show(block=False)
        # plt.show()
        plt.pause(0.1)'''
        plt.close()
        i += 5
        j += 1


def mf_example(waveform: np.array, mf: list, apulse_pos: list):
    shape = mf[0]
    amplitude = mf[1]
    x = [i / 0.64 for i in range(waveform.size)]
    mf_x = [i / 0.64 for i in range(mf[0].size)]

    fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    gs = fig.add_gridspec(2, hspace=0.2)
    axs = gs.subplots()

    axs[0].plot(x, waveform, "k", label="waveform")
    ax_0 = axs[0].twinx()
    ax_0.plot(mf_x, shape, "r", label="mf output")
    ax_0.plot(list(map(mf_x.__getitem__, apulse_pos)),
              list(map(shape.__getitem__, apulse_pos)), "gx")
    ax_0.axhline(0.8, 0, 1800, ls="--", color="g")
    axs[0].legend(loc='lower left')
    ax_0.legend(loc='lower right')
    axs[0].set_xlabel("timestamp /ns")
    axs[0].set_title("MF Shape")
    axs[0].set_ylabel("adc /mV")
    ax_0.set_ylabel("a.u")
    axs[0].grid()
    axs[0].set_xlim(x[256], x[-1])

    axs[1].plot(x, waveform, "k", label="waveform")
    ax_1 = axs[1].twinx()
    ax_1.plot(mf_x, amplitude, "r", label="amplitude")
    ax_1.plot(list(map(mf_x.__getitem__, apulse_pos)),
              list(map(amplitude.__getitem__, apulse_pos)), "gx")
    ax_1.axhline(1, 0, 1800, ls="--", color="g")
    # axs[1].legend(loc='center left')
    ax_1.legend(loc='center right')
    axs[1].set_xlabel("timestamp /ns")
    axs[1].set_title("MF Amplitude")
    axs[1].set_ylabel("adc /mV")
    ax_1.set_ylabel("a.u")
    ax_1.set_ylim(-10, 500)
    axs[1].grid()
    axs[1].set_xlim(x[256], x[-1])

    plt.savefig("plots/mf_example.png")
    plt.close()


def draw_template(om_num: int, om_text: str, template: np.array):
    fig = plt.figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
    x = np.array([i for i in range(template.size)]) / 0.64
    plt.plot(x, template, 'k-')
    plt.xlabel("timestamp /ns")
    plt.ylabel("a.u")
    plt.title(om_text)
    plt.grid()
    plt.savefig(f'plots/template_{om_num}')
    plt.close()


def draw_template_quality(comparison: np.array, templates: list):
    sncalo = calorimeter("template_shape_test")
    sncalo.draw_omid_label()
    sncalo.draw_content_label('{:.3f}')

    for i_temp in range(len(templates)):
        shape = float(np.dot(comparison, templates[i_temp]))
        if shape != 0:
            sncalo.setcontent(i_temp, shape)

    sncalo.draw()
    sncalo.save("plots")


def get_templates(root_file: str, comparison: int):
    template_root_file = ROOT.TFile(root_file, "READ")
    templates = []

    comparison_template = []

    for i in range(712):
        template = []

        h_template = template_root_file.Get(f'Template_Ch{i}')

        try:
            h_template.GetEntries()
        except:
            template = [0 for k in range(templates[0].size)]
            templates.append(np.array(template))
            del h_template
            continue

        for i_bin in range(h_template.GetNbinsX()):
            if i == 1:
                comparison_template.append(h_template.GetBinContent(i_bin))
            template.append(h_template.GetBinContent(i_bin))
        templates.append(np.array(template))
        del h_template

    comparison_template = np.array(comparison_template)
    template_root_file.Close()

    return templates, comparison_template


def draw_AAN(apulse_nums: list):
    sncalo = calorimeter("aan")
    sncalo.draw_omid_label()
    sncalo.draw_content_label('{:.3f}')

    for i in range(len(apulse_nums)):
        aan = float(np.average(apulse_nums[i]))
        if str(aan) == 'nan':
            sncalo.setcontent(i, -1.0)
        else:
            sncalo.setcontent(i, aan)

    sncalo.draw()
    sncalo.save("plots")


def draw_PAR(apulse_nums: list):
    sncalo = calorimeter("par")
    sncalo.draw_omid_label()
    sncalo.draw_content_label('{:.3f}')

    for i in range(len(apulse_nums)):
        if len(apulse_nums[i]) == 0:
            sncalo.setcontent(i, -1)
        else:
            par = float(len(np.where(np.array(apulse_nums[i]) > 0)[0]) / len(apulse_nums[i]) * 100)
            sncalo.setcontent(i, par)

    sncalo.draw()
    sncalo.save("plots")


def draw_ATD(apulse_times: list):
    for i in range(len(apulse_times)):
        canvas = ROOT.TCanvas()
        hist = ROOT.TH1D(f"{om_id_string(i)}", f"{om_id_string(i)}", int((1024 / 0.64) / 20), 0, 1024 / 0.64)
        if np.sum(np.array(apulse_times[i])) == 0:
            continue
        for j in range(len(apulse_times[i])):
            for k in range(apulse_times[i][j]):
                hist.Fill(j/0.64)
        hist.Scale(1.0/hist.GetEntries())
        hist.SetFillColor(2)
        hist.Sumw2()
        hist.SetXTitle("timestamp /ns")
        hist.Draw()
        canvas.SetGrid("HIST")
        canvas.SaveAs(f"plots/h_atd_{i}.png")
        del hist
        del canvas


def draw_AAD(apulse_amplitudes: list):
    for i in range(len(apulse_amplitudes)):
        canvas = ROOT.TCanvas()
        hist = ROOT.TH1D(f"om_id_string(i)", f"{om_id_string(i)}", 40, 0, 400)
        for j in range(len(apulse_amplitudes[i])):
            hist.Fill(apulse_amplitudes[i][j])
        hist.SetXTitle("mf amplitude")
        hist.SetFillColor(2)
        hist.Draw("HIST")
        canvas.SetGrid()
        canvas.SaveAs(f"plots/aad_{i}.png")
        del hist
        del canvas


def draw_event_map(n_events: list):
    sncalo = calorimeter("event_map")
    sncalo.draw_omid_label()
    sncalo.draw_content_label('{}')

    for i in range(len(n_events)):
        sncalo.setcontent(i, n_events[i])

    sncalo.draw()
    sncalo.save("plots")


def main():
    args = io_parse_arguments()
    input_file = args.i
    template_file = args.t

    if input_file is None:
        input_file = "/Users/willquinn/Desktop/test.root"
    if template_file is None:
        template_file = "/Users/willquinn/Desktop/templates.root"

    # templates, comparison_template = get_templates(template_file, comparison=1)
    # draw_template_quality(comparison_template, templates)

    root_file = ROOT.TFile(input_file, "READ")
    tree = root_file.T

    # apulse_nums = [[] for i in range(712)]
    apulse_times = [[0 for i in range(1024)] for j in range(712)]
    # apulse_amplitudes = [[] for i in range(712)]
    # n_events = [0 for i in range(712)]

    i_event = 0
    for event in tree:
        # baseline = event.baseline
        # waveform = np.array(event.waveform)
        OM_ID = event.OM_ID
        # n_events[OM_ID] += 1
        for i_pulse in list(event.apulse_times):
            apulse_times[OM_ID][i_pulse] += 1

        # for i_pulse in list(event.apulse_amplitudes):
        #     apulse_amplitudes[OM_ID].append(i_pulse)
        # apulse_nums[OM_ID].append(event.apulse_num)
        # if event.apulse_num > 4:
        #    sweep_example(om_id_string(OM_ID), templates[OM_ID], waveform-baseline)
        #    break
        # mf_example(waveform-baseline, [np.array(event.mf_shapes), np.array(event.mf_amplitudes)],
        #            list(event.apulse_times))
        i_event += 1
        if not i_event % 100000:
            print(f">>> processed {i_event}/{tree.GetEntries()}")

    # draw_AAN(apulse_nums)
    # draw_PAR(apulse_nums)
    draw_ATD(apulse_times)
    # draw_AAD(apulse_amplitudes)
    # draw_event_map(n_events)

    root_file.Close()


if __name__ == '__main__':
    main()
