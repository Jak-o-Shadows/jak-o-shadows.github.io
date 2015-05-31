# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 13:55:00 2012

@author: Jak_o_Shadows
"""
import sqlite3
import random
import wx


class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.listCtrl = wx.ListCtrl(self, -1, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        self.listCtrl.InsertColumn(0, "a")
        self.listCtrl.InsertColumn(1, "b")
        self.listCtrl.InsertColumn(2, "c")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_LIST_COL_CLICK, self.setupSort, self.listCtrl)
        
        #for allowing sorting in descending order        
        self.oldC = -1
        self.reverse = False
        self.command = "SELECT * FROM tbl"
        
        #populate the database
        self.setupDB()
        
    def __set_properties(self):
        self.SetTitle("Fill a list control with a database!")

    def __do_layout(self):
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()

    def setupSort(self, event):
        """Sets the command for filling the list control, based on
        what column is clicked
        """
        c = event.GetColumn()   #get the column that was clicked on
        
        if c==0:
            #order by first column
            self.command = "SELECT * FROM tbl ORDER BY a"
        elif c==1:
            #order by second column
            self.command = "SELECT * FROM tbl ORDER BY b"
        elif c==2:
            #order by third column
            self.command = "SELECT * FROM tbl ORDER BY c"
                
        
        #Toggle reverse
        if c == self.oldC:
            self.reverse = not self.reverse
        else:
            self.reverse = False
            
        #if reverse, append "DESC" to the select command
        if self.reverse:
            self.command += " DESC"
            
        self.oldC = c
        self.fillLC()
        event.Skip()
        
    def setupDB(self):
        """Open the database, add a table"""
        self.con = sqlite3.connect(":memory:")   #open a connection to a DB in RAM
        command = "CREATE TABLE tbl(id INTEGER PRIMARY KEY AUTOINCREMENT, a INT, b TEXT,\
                    c REAL)"
        self.con.execute(command)    #add a table (tbl) with 4 columns)
        self.con.commit()
        self.fillDB()

    def fillDB(self):
        letters = "abcdefghijklmnopqrstuvwxyz"
        command = "INSERT INTO tbl VALUES(Null, ?, ?, ?)"
        for i in xrange(100):
            a = random.randint(1,26)    #add a random int from 1 to 26
            b = letters[random.randrange(0,26)] #add a random letter from letters
            c = random.random()*26      #add a random float from (0,26)
            self.con.execute(command, (a, b, c))    #add the data to the DB

        #use executemany
        data = []
        for i in xrange(long(1e2)):
            a = random.randint(1,26)    #add a random int from 1 to 26
            b = letters[random.randrange(0,26)] #add a random letter from letters
            c = random.random()*26      #add a random float from (0,26)
            data.append((a, b, c))  #add the data to a list
        self.con.executemany(command, data) #all data in the list is added to db
        self.con.commit()   #commits the data, saving it
        self.fillLC()   #update the list control
        
        
    def fillLC(self):
        """Fills the list control based on the sorting command"""
        self.listCtrl.DeleteAllItems()  #since we're sorting, must delete all
        #then get a list of tuples of all the data
        data = self.con.execute(self.command).fetchall()
        for i in data:
            #loop through and add it
            self.listCtrl.Append(i[1:])   
            
if __name__ == "__main__":
    #wxGlade default stuff
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame_1 = MyFrame(None, -1, "")
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()
