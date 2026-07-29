import time
import serial
from flow import *
from nesp_lib import Port, Pump, PumpingDirection

vrxn = 12*8.75*0.8 # reaction volume under light in ul
vfull = 250 # full cell volume in ul
tbuffer = 30 # time in s after the solution fully fills the cell before recording

tRs = [tR*60 for tR in range(1, 8)] # residence times in s for 1-7 min in 1 min intervals
flowrates = [flowprop(tR=tR/60, vol=vrxn*1e-3) for tR in tRs] # flow rates in ml/min for each tR
points = flowpoints(tRs, tbuffer, vrxn, vfull) # time points in s to record absorbances

print(tRs)
print(flowrates)
print(points)

com = 'COM4'

# pump controls using rs-232 serial commands
def rs232():
    ser = serial.Serial(com, 19200, parity=serial.PARITY_NONE, bytesize=8, stopbits=1, 
                    timeout=1, xonxoff=0, rtscts=0)
    print('Pump connected')

    try:
        def send(cmd):
            ser.write((cmd + '\r\n').encode()) # send commands as encoded strings
            time.sleep(0.05)
            response = ser.readline().decode().strip()
            if response:
                print(f'{cmd} -> {response}')
            else:
                print(f'{cmd} -> (no response)')

        send('RESET')
        send('DIA 9.000')

        for i in range(len(tRs)):
            if i == 0:
                send('PHN 1')
                send('FUN RAT')
                send(f'RAT {flowrates[i]} MM')
                send('DIR INF')
                send('RUN')
            else:
                send(f'RAT {flowrates[i]}')
            time.sleep(points[i])

    except Exception as e:
        print(f'Pump error: {e}')

    finally:
        send('STP')
        print(f'Pump stopped')
        ser.close()

# pump controls using nesp-lib
def nesplib():
    with Port(com) as port:
        pump = Pump(port)
        print('Pump connected')

        try:
            pump.syringe_diameter_mm = 9.000
            pump.pumping_direction = PumpingDirection.INFUSE

            pump.run(False)
            for i in range(len(tRs)):
                pump.pumping_rate_ml_per_min = flowrates[i]
                print(f'New flow rate = {flowrates[i]*1e3:.2f} ul/min for tR = {tRs[i]/60} min')
                time.sleep(points[i])
                
        except Exception as e:
            print(f'Pump error: {e}')

        finally:
            pump.stop()
            print(f'Pump stopped')

rs232()