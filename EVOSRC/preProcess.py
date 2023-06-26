from torch.utils.data import Dataset
import torch
import numpy as np

def addToContext(context, letter):
    newContext = torch.cat((context, letter), dim = 1)[:, letter.shape[1]:]
    return newContext

class chemblData(Dataset):
    def __init__(self, path, maxLen=None, voc=None):
        self.path = path
        fh = open(path, "r")
        self.rawText = fh.read().splitlines()
        fh.close()
        self.smiles = []
        for i in self.rawText:
            self.smiles.append("$" + i.replace("Cl", "D").replace("Br", "E") + "~")
        print(self.smiles[-1])
        if maxLen == None and voc == None:
            self.maxLen, self.voc = self.getMaxLen_voc()
        else:
            self.maxLen = maxLen
            self.voc = voc

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        return self.encode(self.smiles[idx] + "€" * (self.maxLen - len(self.smiles[idx])))

    def getMaxLen_voc(self):
        maxLen = 0
        voc = []
        for i in self.smiles:
            if len(i) > maxLen:
                maxLen = len(i)
            for j in i:
                if j not in voc:
                    voc.append(j)
        voc.append("€")

        vocpath = self.path.split("/")[1][:-4]
        fh = open("data/"+vocpath+".voc", "w")
        fh.write(str(voc))
        fh.close()

        return maxLen, voc

    def encode(self, string):
        outT = torch.zeros(len(string), len(self.voc), dtype=torch.int8)
        for i in range(0, len(string)):
            outT[i, self.voc.index(string[i])] = 1
        return outT


def wordToTensor(word, allLetters):
    """

    Parameters
    ----------
    word : list
        the word in list format(shoul work for other constructs too).
    allLetters : string
        a string containing all the known letters for the problem at hand.

    Returns
    -------
    wordTensor : torch.Tensor
        A len(allLetters) x len(word) tensor containing all letters as one-hot encoded vectors, stacked
        horizontally

    """
    wordTensor = torch.zeros(len(word), len(allLetters))
    for i in range(0, len(word)):
        wordTensor[i] = categoryToTensor(word[i], allLetters).reshape(len(allLetters))
    return wordTensor


def categoryToTensor(category, categoriesList):
    """
    translator in one-hot encoded format

    Parameters
    ----------
    category : whatever (ad esempio una lettera)
        l'oggetto che vuoi trasformare in one-hot encoded format.
    categoriesList : list
        una lista di tutte le categorie possibili per quell'oggetto in modo da
        cotruire l'array.

    Returns
    -------
    cTensor : torch.tensor
        il tensore one-hot. Viene dato come vettore con una sola dimensione
        se sarà riga o colonna poi te la vedi tu

    """
    cTensor = torch.zeros(len(categoriesList))
    cTensor[categoriesList.index(category)] = 1
    return cTensor


def tensorToCategory(tensor, categoriesList):
    """


    Parameters
    ----------
    tensor : torch.tensor
        one-hot encoded tensor.
    categoriesList : list
        list of all possible values for the tensor at hand.

    Returns
    -------
    string (or type of the objects in the list)
        The category label that the input tensor represents.

    """
    index = torch.argmax(tensor)
    return categoriesList[index]


class wordEmbeddedData(Dataset):
    def __init__(self, path, maxLen=None, voc=None):
        self.path = path
        fh = open(path, "r")
        self.rawText = fh.read().splitlines()
        fh.close()
        self.smiles = []
        for i in self.rawText:
            self.smiles.append("$" + i.replace("Cl", "D").replace("Br", "E") + "~")

        if maxLen == None and voc == None:
            self.maxLen, self.voc = self.getMaxLen_voc()
        else:
            self.maxLen = maxLen
            self.voc = voc

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        return self.encode(self.smiles[idx] + "€" * (self.maxLen - len(self.smiles[idx])))

    def getMaxLen_voc(self):
        maxLen = 0
        voc = []
        for i in self.smiles:
            if len(i) > maxLen:
                maxLen = len(i)
            for j in i:
                if j not in voc:
                    voc.append(j)
        voc.append("€")

        vocpath = self.path.split("/")[1][:-4]
        fh = open("data/"+vocpath+".voc", "w")
        fh.write(str(voc))
        fh.close()

        return maxLen, voc

    def encode(self, string):
        outT = torch.zeros(len(string), dtype = torch.float16)
        for i in range(0, len(string)):
            outT[i] = self.voc.index(string[i])
        return outT
