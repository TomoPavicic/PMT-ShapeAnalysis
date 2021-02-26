# include "TBox.h"
# include "TText.h"
# include "TLatex.h"
# include "TLine.h"
# include "TCanvas.h"
# include "TPad.h"
# include "TColor.h"
# include "TSystem.h"
# include "TString.h"

# include<vector>
import ROOT
'''namespace
sndisplay
{'''


class calorimeter:
    '''{
        public:
            calorimeter(const char * n = ""): name(n)
    {'''
    def __init__(self, new_name: str):
        self.name = new_name
        self.canvas_it = None
        self.canvas_fr = None
        #TODO: check this
        self.nb_om = 0
        self.content = []
        self.ombox = []
        self.omid_text_v = []
        self.omnum_text_v = []
        self.content_text_v = []


        self.draw_omid = False
        self.draw_omnum = False
        self.draw_content = False
        self.draw_content_format = "%.0f"

        self.has_italy_data = False
        self.has_french_data = False

        for omnum in range(self.nb_om):
            self.content.append(0)

        self.range_min = -1
        self.range_max = -1

        spacerx = 0.0125
        spacery = 0.0250

        mw_sizey = (1-4 * spacery) / (13+2)
        gv_sizey = mw_sizey
        xw_sizey = mw_sizey * 13. / 16.

        mw_sizex = (1-4 * spacerx) / (20+4)
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
                        x1 += mw_sizex * (19-mw_column)
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
                    omid_text = ROOT.TText(x1+0.5 * mw_sizex, y1+0.667 * mw_sizey, omid_string)
                    omid_text.SetTextSize(0.014)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1+0.5 * mw_sizex, y1+0.333 * mw_sizey, omnum_string)
                    omnum_text.SetTextSize(0.02)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_text = ROOT.TText (x1+0.5 * mw_sizex, y1+0.333 * mw_sizey, "")
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
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (1 - xw_column) * xw_sizex
                            else:
                                x1 = spacerx + xw_sizex * xw_column
                        else:
                            if xw_wall == 0:
                                x1 = spacerx + xw_sizex * xw_column
                            else:
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (1 - xw_column) * xw_sizex

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
                        x1 = spacerx + 2 * xw_sizex + spacerx + gv_sizex * (16-1-gv_column)
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
                    omid_text = ROOT.TText(x1+0.5 * gv_sizex, y1+0.667 * gv_sizey, omid_string)
                    omid_text.SetTextSize(0.014)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1+0.5 * gv_sizex, y1+0.333 * gv_sizey, omnum_string)
                    omnum_text.SetTextSize(0.02)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_text = ROOT.TText (x1+0.5 * gv_sizex, y1+0.333 * gv_sizey, "")
                    content_text.SetTextSize(0.02)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

        it_label = ROOT.TText(spacerx, spacery+gv_sizey+spacery+13 * mw_sizey+spacery+0.25 * gv_sizey, "  ITALY")
        it_label.SetTextSize(0.036)


        fr_label = ROOT.TText (spacerx, spacery+gv_sizey+spacery+13 * mw_sizey+spacery+0.25 * gv_sizey, "FRANCE")
        fr_label.SetTextSize(0.036)

        stops   = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00]
        red     = [0.25, 0.00, 0.20, 1.00, 1.00, 0.90]
        green   = [0.25, 0.80, 1.00, 1.00, 0.80, 0.00]
        blue    = [1.00, 1.00, 0.20, 0.00, 0.00, 0.00]

        nRGBs = 6

        # palette_index = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 100);

    nmwall  = 520
    nxwall  = 128
    ngveto  = 64
    nb_om   = 712

    '''enum
    {
        AUTO,
        FULL,
        FULL_ITALY,
        FULL_FRANCE,
        MW_ITALY,
        MW_FRANCE,
        XW_MOUNTAIN,
        XW_TUNNEL};'''

    def setrange(self, xmin: float, xmax: float):
        self.range_min = xmin
        self.range_max = xmax


    def draw_omid_label(self):
        self.draw_omid = True

    def draw_omnum_label(self):
        self.draw_omnum = True

    def draw_content_label(self, string: str):
        self.draw_content_format = ROOT.TString(string)
        self.draw_content = True

    def draw(self):
        if self.canvas_it == None:
            self.canvas_it = ROOT.TCanvas(f'C_it_{self.name}', self.name, 900, 600)

        if self.canvas_fr == None:
            canvas_fr = ROOT.TCanvas(f'C_fr_{self.name}', name, 900, 600)

        if self.draw_content:
            for omnum in range(self.nb_om):
                ttext = self.content_text_v[omnum]
                ttext.SetText(ttext.GetX(), ttext.GetY(), f'{self.content[omnum]}')


        '''// // // // // // /
        // Draw IT //
        // // // // // // /'''

        self.canvas_it.cd()
        self.canvas_it.SetEditable(True)

        mw_side=0
        for mw_column in range(20):
            for mw_row in range(13):
                id = mw_side * 20 * 13 + mw_column * 13 + mw_row
                self.ombox[id].Draw("l")
                if (self.draw_omid):
                    omid_text_v[id]->Draw();
                if (draw_omnum) omnum_text_v[id]->Draw();
                if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
                // else if (draw_omid_num) omid_num_text_v[id]->Draw();
            }
        }

    unsigned int xw_side=0;
    for (unsigned int xw_wall=0; xw_wall < 2; ++xw_wall) {
    for (unsigned int xw_column=0; xw_column < 2; ++xw_column) {
    for (unsigned int xw_row=0; xw_row < 16; ++xw_row) {
    unsigned int id = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row;
    ombox[id]->Draw("l");
    if (draw_omid) // if ((xw_column % 5) == 0)
    omid_text_v[id]->Draw();
    if (draw_omnum) omnum_text_v[id]->Draw();
    if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
    // else if (draw_omid_num) omid_num_text_v[id]->Draw();
    }
    }
    }

    unsigned int gv_side=0;
    for (unsigned int gv_wall=0; gv_wall < 2; ++gv_wall) {
    for (unsigned int gv_column=0; gv_column < 16; ++gv_column) {
    unsigned int id = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column;
    ombox[id]->Draw("l");
    if (draw_omid) // if ((gv_column % 5) == 0)
    omid_text_v[id]->Draw();
    if (draw_omnum) omnum_text_v[id]->Draw();
    if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
    // else if (draw_omid_num) omid_num_text_v[id]->Draw();
    }
    }

    it_label->Draw();

    canvas_it->SetEditable(False);

    // // // // // // /
    // Draw FR //
    // // // // // // /

    canvas_fr->cd();
    canvas_fr->SetEditable(True);

    mw_side=1;
    for (unsigned int mw_column=0; mw_column < 20; ++mw_column) {
    for (unsigned int mw_row=0; mw_row < 13; ++mw_row) {
    unsigned int id = mw_side * 20 * 13 + mw_column * 13 + mw_row;
    ombox[id]->Draw("l");
    if (draw_omid) // if (((mw_column % 5) == 0) | | (mw_column == 19))
    omid_text_v[id]->Draw();
    if (draw_omnum) omnum_text_v[id]->Draw();
    if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
    // else if (draw_omid_num) omid_num_text_v[id]->Draw();
    }
    }

    xw_side=1;
    for (unsigned int xw_wall=0; xw_wall < 2; ++xw_wall) {
    for (unsigned int xw_column=0; xw_column < 2; ++xw_column) {
    for (unsigned int xw_row=0; xw_row < 16; ++xw_row) {
    unsigned int id = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row;
    ombox[id]->Draw("l");
    if (draw_omid) // if ((xw_column % 5) == 0)
    omid_text_v[id]->Draw();
    if (draw_omnum) omnum_text_v[id]->Draw();
    if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
    // else if (draw_omid_num) omid_num_text_v[id]->Draw();
    }
    }
    }

    gv_side=1;
    for (unsigned int gv_wall=0; gv_wall < 2; ++gv_wall) {
    for (unsigned int gv_column=0; gv_column < 16; ++gv_column) {
    unsigned int id = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column;
    ombox[id]->Draw("l");
    if (draw_omid) // if ((gv_column % 5) == 0)
    omid_text_v[id]->Draw();
    if (draw_omnum) omnum_text_v[id]->Draw();
    if (draw_content & & content[id] != 0) content_text_v[id]->Draw();
    // else if (draw_omid_num) omid_num_text_v[id]->Draw();
    }
    }

    fr_label->Draw();
    canvas_fr->SetEditable(False);

    //

    update();
    }


void reset() {
for (unsigned int omnum=0; omnum < nb_om; ++omnum)
content[omnum] = 0;

// float content_min = 0; float content_max = 1;

for (unsigned int omnum=0; omnum < nb_om; ++omnum)
ombox[omnum]->SetFillColor(0);

canvas_it->Modified();
canvas_it->Update();
gSystem->ProcessEvents();

}

float getcontent (int omnum)
{
return content[omnum];
}

void
setcontent(int
omnum, float
value)
{
if (omnum >= 0 & & omnum < 260)
    has_italy_data = True;
else if (omnum < 520)
has_french_data = True;
else if (omnum < 584)
    has_italy_data = True;
else if (omnum < 648)
has_french_data = True;
else if (omnum < 680)
    has_italy_data = True;
else if (omnum < 712)
has_french_data = True;

content[omnum] = value;
}

void
setcontent(int
om_side, int
om_wall, int
om_column, int
om_row, float
value)
{
int
omnum = -1;

// auto
detect
MW
if ((om_side != -1) & & (om_wall == -1) & & (om_column != -1) & & (om_row != -1))
    omnum = 260 * om_side + 13 * om_column + om_row;

// auto
detect
XW
else if ((om_side != -1) & & (om_wall != -1) & & (om_column != -1) & & (om_row != -1))
    omnum = 520 + 64 * om_side + 32 * om_wall + 16 * om_column + om_row;

// auto
detect
GV
else if ((om_side != -1) & & (om_wall != -1) & & (om_column != -1) & & (om_row == -1))
    omnum = 520 + 128 + 32 * om_side + 16 * om_wall + om_column;

else {
printf("+++ sndisplay: skipping OM (%d.%d.%d.%d)\n", om_side, om_wall, om_column, om_row);
return;}

if (om_side == 0)
has_italy_data = True;
else if (om_side == 1)
has_french_data = True;

content[omnum] = value;
}


void fill (int omnum, float value=1)
{
setcontent(omnum, content[omnum]+value);
}

// void fill (int om_side, int om_wall, int om_column, int om_row, float value=1)
// {
// setcontent(om_side, om_wall, om_column, om_row, content[omnum]+value);
//}

void update()
{
float content_min = content[0];
float content_max = content[0];

if (range_min == -1) {
for (unsigned int omnum=1; omnum < nb_om; ++omnum) {
if (content[omnum] < content_min) content_min = content[omnum];
if (content[omnum] > content_max) content_max = content[omnum];}
} else {
range_min = 0;
content_max = range_max;}

for (unsigned int omnum=0; omnum < nb_om; ++omnum)
if (content[omnum] != 0)
   // ombox[omnum]->SetFillColor(TColor::
    GetColorPalette((int)(99 * (content[omnum] - content_min) / (content_max - content_min))));
ombox[omnum]->SetFillColor(palette_index + (int)(99 * (content[omnum] - content_min) / (content_max - content_min)));
else
ombox[omnum]->SetFillColor(0);

canvas_it->Modified();
canvas_it->Update();

canvas_fr->Modified();
canvas_fr->Update();

gSystem->ProcessEvents();
}

// void
setcontent(geomid)

// void
setcontent(uint32_t
geomid0, uint32_t
geomid1, uint32_t
geomid2, uint32_t
geomid3, uint32_t
geomid4, float
content)
// {
   // switch(geomid0)
   // {
// case
// }
// }

TString name;

float range_min, range_max;

std::vecto r <floa t> content;
std::vecto r <TBo x *>   ombox;

bool draw_omid;
bool draw_omnum;
bool draw_content;
TString draw_content_format;

bool has_italy_data;
bool has_french_data;

std::vecto r <TTex t *> omid_text_v;
std::vecto r <TTex t *> omnum_text_v;
std::vecto r <TTex t *> content_text_v;

TLine *source_foil;
TText * it_label;
TText * fr_label;

// TCanvas * canvas;
TCanvas * canvas_it;
TCanvas * canvas_fr;

TPad * pad_italy;
TPad * pad_french;

int
palette_index;
};
}

void
sndisplay_test()
{
    sndisplay:: calorimeter * sncalo = new
sndisplay::calorimeter;

sncalo->draw_omid_label();
sncalo->draw_content_label("%.1f");

TRandom
trand;

for (int omnum=0; omnum < 520; ++omnum)
sncalo->setcontent(omnum, trand.Gaus(100, 10));

for (int omnum=520; omnum < 712; ++omnum)
sncalo->setcontent(omnum, trand.Gaus(50, 10));

sncalo->draw();

}