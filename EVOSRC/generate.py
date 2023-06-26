import torch
import net
import preProcess
from preProcess import chemblData
import torch.nn as nn
import random
import sys
import rdkit.Chem as Chem
import rdkit.Chem.rdchem
import SAScore
import os


from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


class generatorWithSubstitutions():
    def __init__(self, paramsDict):
        self.device = paramsDict["device"]
        self.name = paramsDict["name"]
        self.network = torch.load(paramsDict["networkPath"], map_location=torch.device('cpu')).to(self.device)
        fh = open(paramsDict["vocabularyPath"])
        self.voc = eval(fh.read())
        fh.close()
        self.contextLength= paramsDict["contextLength"]
        self.vocLength = len(self.voc)
        self.softmax = nn.Softmax(dim = -1)
        self.rawInputMol = paramsDict["inputString"]
        self.inputMole = "$" + paramsDict["inputString"].replace("Cl", "D").replace("Br", "E") + "~"
        self.wantedValids = paramsDict["requestedValid"]
        #Qesto parametro va tenuto d'occhio!
        self.N = self.wantedValids*1000
        self.queryMol = rdkit.Chem.MolFromSmiles(self.rawInputMol)
        self.querySA = SAScore.calculateScore(self.queryMol)
        self.positions = params["positions"]

    def run(self):
        thisPath = os.getcwd()
        query = preProcess.wordToTensor(self.inputMole, self.voc)
        qs = query.shape[0]
        vocLen = len(self.voc)
        for z in range(0, qs):
            if preProcess.tensorToCategory(query[z], self.voc) == "~":
                maxQ = z
                break
        molFile = ""
        j = 0
        valids = 0
        while j < self.N and valids < self.wantedValids:
            outT = 0
            context = torch.zeros(1,self.contextLength*len(self.voc)).to(self.device)
            cellT1 = torch.zeros(1, self.network.n_h1, device = self.device)
            cellT2 = torch.zeros(1, self.network.n_h2, device = self.device)
            hiddenT1 = torch.zeros(1, self.network.n_h1, device = self.device)
            hiddenT2 = torch.zeros(1, self.network.n_h2, device = self.device)
            answer = "$"
            context = preProcess.addToContext(context, preProcess.categoryToTensor("$", self.voc).unsqueeze(0).to(self.device))
        #   #Determino le posizioni da sostituire
        # Not today!

        
            for letter in range(1, qs):
                if letter in self.positions:
                    charT, hiddenT1, hiddenT2, cellT1, cellT2 = self.network(context, hiddenT1, hiddenT2, cellT1, cellT2)
                    charTemp = self.softmax(charT)
                    s = 0
                    rn = random.uniform(0,1)
                    for vocPosition in range(0, vocLen):
                        maxv, maxp = torch.max(charTemp, dim = -1)
                        s = s + maxv.item()
                        if s >= rn:
                            p = maxp.item()
                            break
                        else:
                            charTemp[0][maxp.item()] = 0

                    outT = torch.zeros(len(self.voc)).to(self.device)
                    outT[p] = 1
                    context = preProcess.addToContext(context, outT.unsqueeze(0))

                    answer = answer + preProcess.tensorToCategory(outT, self.voc)

                    if preProcess.tensorToCategory(outT, self.voc) == "~" or len(answer)>=200:
                        break
                else:
                    answer = answer + preProcess.tensorToCategory(query[letter], self.voc)
                    charT, hiddenT1, hiddenT2, cellT1, cellT2 = self.network(context, hiddenT1, hiddenT2, cellT1, cellT2)
                    context = preProcess.addToContext(context, query[letter].unsqueeze(0))

                    if preProcess.tensorToCategory(query[letter], self.voc) == "~" or len(answer)>=200:
                        break

            answer = answer.replace("D", "Cl").replace("E", "Br").replace("€", "").replace("~", ""). replace("$", "")
    
            j = j + 1
            
            try:
                Mol = 0
                mol = Chem.MolFromSmiles(answer)
                valid = True
                if mol == None:
                    valid = False
                Mol = rdkit.Chem.rdchem.Mol(mol)
                rings = Mol.GetRingInfo().AtomRings()
                for i in range(0, len(rings)):
                    if len(rings[i]) > 8:
                        valid = False
                        break
                try:
                    SA = SAScore.calculateScore(mol)
                except:
                    SA = 10
                if self.querySA >= 3:
                    threshhold = self.querySA+1
                else:
                    threshhold = 4

                if SA>=threshhold:
                    valid =False

            except Exception as e:
                valid = False
            
            if valid ==True:
                answer1 = Chem.MolToSmiles(mol)
                if answer1 not in molFile and answer1!=self.rawInputMol:
                    valids += 1
                    molFile += answer1 + "\n"

        molFile = molFile.replace("D", "Cl").replace("E", "Br").replace("€", "").replace("~", ""). replace("$", "")

        
        
        fh = open(thisPath + "/" + self.name + ".csv", "w")
        fh.write(molFile)
        fh.close()
        
        return molFile


scriptPath = os.path.dirname(os.path.realpath(__file__)) 

params={"name": sys.argv[2],
         "vocabularyPath": scriptPath + "/data/CHEMBL28-0.05.voc",
         "networkPath": scriptPath + "/trainedNetworks/trained2021-11-22.pt",
         "contextLength": 8,
         "device" : "cpu",
         "inputString": sys.argv[1],
         "requestedValid": int(sys.argv[3]),
         "positions": eval(sys.argv[4])
         }

gen = generatorWithSubstitutions(params)
print(gen.run())
