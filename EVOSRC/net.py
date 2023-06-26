"""
This module contains the building blocks for my nets.
"""


import torch.nn as nn
import torch
import preProcess
import numpy as np

class myLSTM(nn.Module):
    def __init__(self, n_i, n_h):
        super(myLSTM, self).__init__()
        self.n_h = n_h
        self.n_i = n_i

        self.f = nn.Linear(self.n_i + self.n_h, self.n_h)
        torch.nn.init.xavier_uniform_(self.f.weight, gain = 1)
        self.sigmoid = nn.Sigmoid()

        self.i = nn.Linear(self.n_i + self.n_h, self.n_h)
        torch.nn.init.xavier_uniform_(self.i.weight, gain = 1)

        self.c = nn.Linear(self.n_i + self.n_h, self.n_h)
        torch.nn.init.xavier_uniform_(self.c.weight, gain = 5/3)

        self.tanh = nn.Tanh()


        self.o = nn.Linear(self.n_i + self.n_h, self.n_h)
        torch.nn.init.xavier_uniform_(self.o.weight, gain = 1)


    def forward(self, X, hidden, cell):

        concatInput = torch.cat((X,hidden), 1)

        f = self.sigmoid(self.f(concatInput))
        i = self.sigmoid(self.i(concatInput))
        c = self.tanh(self.c(concatInput))
        o = self.sigmoid(self.o(concatInput))
        h = o * self.tanh(c)
        out = h

        cellT = f * cell + c * i

        return out, h, cellT




class QueST(nn.Module):
    def __init__(self, n_i, n_h1, n_h2, n_h3, n_h4):
        super(QueST, self).__init__()
        self.n_h1 = n_h1
        self.n_i = n_i
        self.n_h2 = n_h2

        self.lstm1 = myLSTM(n_i, n_h1)
        self.lstm2 = myLSTM(n_h1, n_h2)
        self.linear3 = nn.Linear(n_h2,n_h3)
        torch.nn.init.xavier_uniform_(self.linear3.weight, gain = np.sqrt(2))
        self.relu = nn.ReLU()
        self.linear4 =nn.Linear(n_h3, n_h4)
        torch.nn.init.xavier_uniform_(self.linear4.weight, gain = 1)

    def forward(self, X, hidden1, hidden2, cell1, cell2):
        out1 , h1, c1 = self.lstm1(X, hidden1, cell1)
        out2 ,h2, c2 = self.lstm2(out1, hidden2, cell2)
        out3 = self.relu(self.linear3(out2))
        out = self.linear4(out3)

        return out, h1, h2,  c1, c2



# net = QueST(37, 512, 512, 256, 37)
# testInp = torch.zeros(64,37)
# testInp[:, 5] = 1
# h1 = torch.zeros(64,512)
# h2 = torch.zeros(64,512)
# c1 = torch.zeros(64,512)
# c2 = torch.zeros(64,512)
# out, h1, h2, c1, c2 = net(testInp, h1, h2, c1, c1)
# print(out.shape)

# =============================================================================
# fh = open("/media/peppe/e9674da5-085c-4d33-a303-409db332053c/Cose/PhD/traduzionePyTorch/data/vocabulary")
# voc = eval(fh.read())
# fh.close()
# =============================================================================
