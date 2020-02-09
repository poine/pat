#! /usr/bin/env python
import sys, os, time, math, numpy as np
import logging
import re
import pdb
import matplotlib.pyplot as plt
import pprzlink.message
import pprzlink.messages_xml_map


def parse_msg(line):
    return re.findall('^(\S*) (\S*) (\S*) (.*)$', line)[0] # ts, sender, msg_id, payload

class PprzBag:
    def __init__(self, filename):
        self.filename = filename

    def read_messages(self, msg_class, msg_ids):
        with open(self.filename) as fp:
            msgs = {msg_id:pprzlink.message.PprzMessage(msg_class, msg_id) for msg_id in msg_ids}
            for cnt, line in enumerate(fp):
                ts, sender, msg_id, payload = parse_msg(line)
                try:
                    msgs[msg_id].ivy_string_to_payload(payload)
                    yield float(ts), int(sender), msgs[msg_id]
                except KeyError: pass
                

def test(filename):
    bag = PprzBag(filename)
    msg_class, msg_ids = 'telemetry', ['AIRSPEED_MS45XX']
    for ts, sender, msg in bag.read_messages(msg_class, msg_ids): 
        print ts, sender, msg

def plot_airspeed(filename):
    time, airspeeds = [],[]
    msg_class, msg_ids = 'telemetry', ['AIRSPEED_MS45XX']
    for ts, sender, msg in PprzBag(filename).read_messages(msg_class, msg_ids): 
        time.append(ts); airspeeds.append(float(msg.airspeed))
    print('read {} messages'.format(len(time)))
    time, airspeeds = np.asarray(time), np.asarray(airspeeds)
    #which = slice(0,20000,1)
    which = time>200
    plt.plot(time[which], airspeeds[which])   


class DataSet:
    def __init__(self):
        self.msg_class, self.msg_ids = 'telemetry', ['GPS', 'ACTUATORS', 'IMU_GYRO']
        def rec_gps(msg): return [float(msg.alt)/1000, float(msg.climb)/100]
        def rec_actuators(msg): return [float(_v) for _v in msg.values]
        def rec_imu_gyro(msg): return [float(_v) for _v in [msg.gp, msg.gq, msg.gr]]
        self.cbks = [rec_gps, rec_actuators, rec_imu_gyro]
        
    def readPprzBag(self, filename):
        self.tss = [[] for _m in self.msg_ids]
        self.vals = [[] for _m in self.msg_ids]
        msg_id_of_name = {msg_id:_i for _i, msg_id in enumerate(self.msg_ids)}
        for ts, sender, msg in PprzBag(filename).read_messages(self.msg_class, self.msg_ids):
            _id = msg_id_of_name[msg.name]
            self.tss[_id].append(ts)
            self.vals[_id].append(self.cbks[_id](msg))
    
    def save(self, filename):
        print('saving to {}'.format(filename))
        np.savez(filename, tss=np.asarray(self.tss), vals=np.asarray(self.vals))

    def load(self, filename):
        print('loading from {}'.format(filename))
        data =  np.load(filename, allow_pickle=True)
        self.tss, self.vals = data['tss'], data['vals']

    

def fit_actuators(filename):
    d = DataSet()
    start_time = time.time()
    if 1:
        d.readPprzBag(filename)
        d.save('/tmp/foo.npz')
    else:
        d.load('/tmp/foo.npz')
    elapsed_time = time.time() - start_time
    print elapsed_time
    plt.plot(d.tss[0], d.vals[0])
    plt.figure()
    plt.subplot(2,1,1)
    pdb.set_trace()
    plt.plot(d.tss[1], d.vals[1][:,0])
    plt.subplot(2,1,1)
    plt.plot(d.tss[1], d.vals[1][:,1])
    return
   
    
    
def main(filename):
    #test(filename)
    plot_airspeed(filename)
    #fit_actuators(filename)
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    #filename = '/home/poine/work/FlightEx/Flight_Data/18_06_01__10_51_00_SD.data'
    filename = '/home/poine/work/FlightEx/Flight_Data/19_03_07__11_22_22_SD.data'
    main(filename)
    
