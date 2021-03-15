import ROOT
import numpy as np
from random import gauss


class calorimeter:

    def __init__(self, new_name: str):
        self.name = new_name
        self.canvas_it = ROOT.TCanvas(f'C_it_{self.name}', self.name, 900, 600)
        self.canvas_fr = ROOT.TCanvas(f'C_fr_{self.name}', self.name, 900, 600)
        # TODO: check this
        self.nmwall = 520
        self.nxwall = 128
        self.ngveto = 64
        self.nb_om = 712
        self.content = []
        self.ombox = []
        self.omid_text_v = []
        self.omnum_text_v = []
        self.content_text_v = []

        self.draw_omid = False
        self.draw_omnum = False
        self.draw_content = False
        self.draw_content_format = '{:.0f}'

        self.has_italy_data = False
        self.has_french_data = False

        for omnum in range(self.nb_om):
            self.content.append(None)

        self.range_min = -1
        self.range_max = -1

        spacerx = 0.0125
        spacery = 0.0250

        mw_sizey = (1 - 4 * spacery) / (13 + 2)
        gv_sizey = mw_sizey
        xw_sizey = mw_sizey * 13. / 16.

        mw_sizex = (1 - 4 * spacerx) / (20 + 4)
        gv_sizex = mw_sizex * 20. / 16.
        xw_sizex = mw_sizex

        '''// // // // // // // // // // // // //
        // MWALL initialisation //
        // // // // // // // // // // // // //'''

        for mw_side in range(2):

            for mw_column in range(20):

                for mw_row in range(13):
                    omnum = mw_side * 20 * 13 + mw_column * 13 + mw_row

                    x1 = spacerx + 2 * xw_sizex + spacerx

                    if mw_side == 0:
                        x1 += mw_sizex * (19 - mw_column)
                    else:
                        x1 += mw_sizex * mw_column

                    y1 = spacery + gv_sizey + spacery + mw_sizey * mw_row
                    x2 = x1 + mw_sizex
                    y2 = y1 + mw_sizey

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.ombox.append(box)

                    omid_string = f'M:{mw_side}.{mw_column}.{mw_row}'
                    omid_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.667 * mw_sizey, omid_string)
                    omid_text.SetTextSize(0.014)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.333 * mw_sizey, omnum_string)
                    omnum_text.SetTextSize(0.02)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.333 * mw_sizey, "")
                    content_text.SetTextSize(0.02)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

        '''// // // // // // // // // // // // //
        // XWALL initialisation //
        // // // // // // // // // // // // //'''

        for xw_side in range(2):

            for xw_wall in range(2):

                for xw_column in range(2):

                    for xw_row in range(16):

                        omnum = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row

                        if xw_side == 0:
                            if xw_wall == 0:
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (
                                            1 - xw_column) * xw_sizex
                            else:
                                x1 = spacerx + xw_sizex * xw_column
                        else:
                            if xw_wall == 0:
                                x1 = spacerx + xw_sizex * xw_column
                            else:
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (
                                            1 - xw_column) * xw_sizex

                        x2 = x1 + xw_sizex
                        y1 = spacery + gv_sizey + spacery + xw_sizey * xw_row
                        y2 = spacery + gv_sizey + spacery + xw_sizey * (xw_row + 1)

                        box = ROOT.TBox(x1, y1, x2, y2)
                        box.SetFillColor(0)
                        box.SetLineWidth(1)
                        self.ombox.append(box)

                        omid_string = f'X:{xw_side}.{xw_wall}.{xw_column}.{xw_row}'
                        omid_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.6 * mw_sizey, omid_string)
                        omid_text.SetTextSize(0.014)
                        omid_text.SetTextAlign(22)
                        self.omid_text_v.append(omid_text)

                        omnum_string = f'{omnum}'
                        omnum_text = ROOT.TText(x1 + 0.5 * xw_sizex, y1 + 0.333 * xw_sizey, omnum_string)
                        omnum_text.SetTextSize(0.02)
                        omnum_text.SetTextAlign(22)
                        self.omnum_text_v.append(omnum_text)

                        content_text = ROOT.TText(x1 + 0.5 * xw_sizex, y1 + 0.333 * xw_sizey, "")
                        content_text.SetTextSize(0.02)
                        content_text.SetTextAlign(22)
                        self.content_text_v.append(content_text)

        '''// // // // // // // // // // // // //
        // GVETO initialisation //
        // // // // // // // // // // // // //'''

        for gv_side in range(2):

            for gv_wall in range(2):

                for gv_column in range(16):

                    omnum = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column

                    if gv_side == 0:
                        x1 = spacerx + 2 * xw_sizex + spacerx + gv_sizex * (16 - 1 - gv_column)
                    else:

                        x1 = spacerx + 2 * xw_sizex + spacerx + gv_sizex * gv_column

                    x2 = x1 + gv_sizex
                    y1 = spacery + gv_wall * (gv_sizey + spacery + 13 * mw_sizey + spacery)
                    y2 = y1 + gv_sizey

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.ombox.append(box)

                    omid_string = f'G:{gv_side}.{gv_wall}.{gv_column}'
                    omid_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.667 * gv_sizey, omid_string)
                    omid_text.SetTextSize(0.014)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.333 * gv_sizey, omnum_string)
                    omnum_text.SetTextSize(0.02)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.333 * gv_sizey, "")
                    content_text.SetTextSize(0.02)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

        self.it_label = ROOT.TText(spacerx, spacery + gv_sizey + spacery + 13 * mw_sizey + spacery + 0.25 * gv_sizey,
                                   "  ITALY")
        self.it_label.SetTextSize(0.036)

        self.fr_label = ROOT.TText(spacerx, spacery + gv_sizey + spacery + 13 * mw_sizey + spacery + 0.25 * gv_sizey,
                                   "FRANCE")
        self.fr_label.SetTextSize(0.036)

        stops = np.array([0.00, 0.20, 0.40, 0.60, 0.80, 1.00], dtype='float64')
        red = np.array([0.25, 0.00, 0.20, 1.00, 1.00, 0.90], dtype='float64')
        green = np.array([0.25, 0.80, 1.00, 1.00, 0.80, 0.00], dtype='float64')
        blue = np.array([1.00, 1.00, 0.20, 0.00, 0.00, 0.00], dtype='float64')

        nRGBs = 6

        self.palette_index = ROOT.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 100)

    '''enum
    {
        AUTO,
        FULL,
        FULL_ITALY,
        FULL_FRANCE,
        MW_ITALY,
        MW_FRANCE,
        XW_MOUNTAIN,
        XW_TUNNEL}'''

    def setrange(self, xmin: float, xmax: float):
        self.range_min = xmin
        self.range_max = xmax

    def draw_omid_label(self):
        self.draw_omid = True

    def draw_omnum_label(self):
        self.draw_omnum = True

    def draw_content_label(self, string: str):
        self.draw_content_format = string
        self.draw_content = True

    def draw(self):

        if self.draw_content:
            for omnum in range(self.nb_om):
                if self.content[omnum] is not None:
                    ttext = self.content_text_v[omnum]
                    ttext.SetText(ttext.GetX(), ttext.GetY(), self.draw_content_format.format(self.content[omnum]))

        '''// // // // // // /
        // Draw IT //
        // // // // // // /'''

        self.canvas_it.cd()
        self.canvas_it.SetEditable(True)

        mw_side = 0
        for mw_column in range(20):
            for mw_row in range(13):
                id = mw_side * 20 * 13 + mw_column * 13 + mw_row
                self.ombox[id].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id].Draw()
                if self.draw_content and self.content[id] is not None:
                    self.content_text_v[id].Draw()

        xw_side = 0
        for xw_wall in range(2):
            for xw_column in range(2):
                for xw_row in range(16):
                    id = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row
                    self.ombox[id].Draw("l")
                    if self.draw_omid:
                        self.omid_text_v[id].Draw()
                    if self.draw_omnum:
                        self.omnum_text_v[id].Draw()
                    if self.draw_content and self.content[id] is not None:
                        self.content_text_v[id].Draw()

        gv_side = 0
        for gv_wall in range(2):
            for gv_column in range(16):
                id = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column
                self.ombox[id].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id].Draw()
                if self.draw_content and self.content[id] is not None:
                    self.content_text_v[id].Draw()

        self.it_label.Draw()
        self.canvas_it.SetEditable(False)

        '''// // // // // // /
        // Draw FR //
        // // // // // // /'''

        self.canvas_fr.cd()
        self.canvas_fr.SetEditable(True)

        mw_side = 1
        for mw_column in range(20):
            for mw_row in range(13):
                id = mw_side * 20 * 13 + mw_column * 13 + mw_row
                self.ombox[id].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id].Draw()
                if self.draw_content and self.content[id] is not None:
                    self.content_text_v[id].Draw()

        xw_side = 1
        for xw_wall in range(2):
            for xw_column in range(2):
                for xw_row in range(16):
                    id = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row
                    self.ombox[id].Draw("l")
                    if self.draw_omid:
                        self.omid_text_v[id].Draw()
                    if self.draw_omnum:
                        self.omnum_text_v[id].Draw()
                    if self.draw_content and self.content[id] is not None:
                        self.content_text_v[id].Draw()

        gv_side = 1
        for gv_wall in range(2):
            for gv_column in range(16):
                id = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column
                self.ombox[id].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id].Draw()
                if self.draw_content and self.content[id] is not None:
                    self.content_text_v[id].Draw()

        self.fr_label.Draw()
        self.canvas_fr.SetEditable(False)

        self.update()

    def reset(self):
        for omnum in range(self.nb_om):
            self.content[omnum] = None

        for omnum in range(self.nb_om):
            self.ombox[omnum].SetFillColor(0)

        self.canvas_it.Modified()
        self.canvas_it.Update()
        ROOT.gSystem.ProcessEvents()

    def getcontent(self, omnum: int):
        return self.content[omnum]

    def setcontent(self, omnum: int, value: float):
        if 0 <= omnum < 260:
            self.has_italy_data = True
        elif omnum < 520:
            self.has_french_data = True
        elif omnum < 584:
            self.has_italy_data = True
        elif omnum < 648:
            self.has_french_data = True
        elif omnum < 680:
            self.has_italy_data = True
        elif omnum < 712:
            self.has_french_data = True

        self.content[omnum] = value

    def setcontent_(self, om_side: int, om_wall: int, om_column: int, om_row: int, value: float):
        omnum = -1

        # auto detect MW
        if om_side != -1 and om_wall == -1 and om_column != -1 and om_row != -1:
            omnum = 260 * om_side + 13 * om_column + om_row
        # auto detect XW
        elif om_side != -1 and om_wall != -1 and om_column != -1 and om_row != -1:
            omnum = 520 + 64 * om_side + 32 * om_wall + 16 * om_column + om_row
        # auto detect GV
        elif om_side != -1 and om_wall != -1 and om_column != -1 and om_row == -1:
            omnum = 520 + 128 + 32 * om_side + 16 * om_wall + om_column
        else:
            print(f'+++ sndisplay: skipping OM ({om_side}.{om_wall}.{om_column}.{om_row})\n')
            return

        if om_side == 0:
            self.has_italy_data = True
        elif om_side == 1:
            self.has_french_data = True

        self.content[omnum] = value

    def fill(self, omnum: int):
        if self.content[omnum] is None:
            val = 1
        else:
            val = self.content[omnum] + 1
        self.setcontent(omnum, val)

    def update(self):

        for i_om in range(self.nb_om):
            if self.content[i_om] is not None:
                content_min = self.content[i_om]
                content_max = self.content[i_om]
                break

        if self.range_min == -1:
            for omnum in range(1, self.nb_om):
                if self.content[omnum] is None:
                    continue
                if self.content[omnum] < content_min:
                    content_min = self.content[omnum]
                if self.content[omnum] > content_max:
                    content_max = self.content[omnum]
        else:
            self.range_min = 0
            content_max = self.range_max

        for omnum in range(self.nb_om):
            if self.content[omnum] is not None:
                self.ombox[omnum].SetFillColor(
                    self.palette_index + int(99 * (self.content[omnum] - content_min) / (content_max - content_min)))
            else:
                self.ombox[omnum].SetFillColor(14)

        self.canvas_it.Modified()
        self.canvas_it.Update()

        self.canvas_fr.Modified()
        self.canvas_fr.Update()

        ROOT.gSystem.ProcessEvents()

    def save(self, location: str):
        self.canvas_it.SaveAs(location + "/" + self.name + '_it.png')
        self.canvas_fr.SaveAs(location + "/" + self.name + '_fr.png')


def sndisplay_test():
    sncalo = calorimeter('test')

    sncalo.draw_omid_label()
    sncalo.draw_content_label('{:.1f}')

    for omnum in range(520):
        sncalo.setcontent(omnum, gauss(100.0, 10.0))

    '''for omnum in range(712):
        sncalo.setcontent(omnum, gauss(50.0, 10.0))'''

    sncalo.draw()
    sncalo.save(".")


if __name__ == '__main__':
    sndisplay_test()
