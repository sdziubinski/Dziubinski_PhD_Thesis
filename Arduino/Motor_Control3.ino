#include <Stepper.h>
#include <string.h>

//constants
const int stepsPerRevolution = 200;
const double conversion = 0.025;

//actual distance between PMTs
const double PMT_separation = 10.0;
//distance between steps when scanning the PMT surface
//const double PMT_separation = 2.0; 

//short boards
const double first_position = 19.6 - 1.84;
//long boards
//const double first_position = 20.6;

//variables
double position = 0.0;
int steps = 0;
int docked = 0;

//pins
const int extreme_LS = 2;
const int home_LS = 3;
const int BACK = 5;
const int power = 4;
const int ENB = 6;
const int FWD = 7;
const int IN1 = 8;
const int IN2 = 9;
const int IN3 = 10;
const int IN4 = 11;
const int RUN = 13;
const int HOME = 12;

Stepper myStepper(stepsPerRevolution, 8, 9, 10, 11);


void setup() {
  //speed control (need to tune to load)
  myStepper.setSpeed(100);

  //initialize the serial port
  Serial.begin(9600);
  Serial.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
  Serial.println("Code Starting...");

  //initialize pins
  pinMode(home_LS, INPUT);
  pinMode(extreme_LS, INPUT);
  pinMode(RUN, INPUT);
  pinMode(HOME, INPUT);
  pinMode(BACK, INPUT);
  pinMode(FWD, INPUT);
  pinMode(power, OUTPUT);

  attachInterrupt(digitalPinToInterrupt(home_LS), home, FALLING);
  //attachInterrupt(digitalPinToInterrupt(extreme_LS), extreme, FALLING);

  digitalWrite(power, HIGH);
}

void loop() {
  // put your main code here, to run repeatedly:
  delay(200);
  //Serial.println(digitalRead(extreme_LS));
  docked = 0;
  /*if (digitalRead(home_LS)==0){
    position = -20.6;
  }
  else {
    docked=0;
  }
    */
  //Serial.println(digitalRead(home_LS));
  if (digitalRead(BACK)) {
    //Serial.print("Moving backward from "); //PMT Gain Correction Stand
    Serial.print("Moving up ");
    //Serial.print(position);
    //Serial.println(" mm.");
    move_backward(PMT_separation);
  }
  else if (digitalRead(FWD)) {
    Serial.print("Moving forward from ");
    Serial.print(position);
    Serial.println(" mm.");
    move_forward(PMT_separation);
  }
  else if(digitalRead(HOME)) {
    Serial.println(digitalRead(home_LS));
    //am I home already?
    if (digitalRead(home_LS) == 0) {
      Serial.println("Already Home");
    }
    else {
      Serial.println("Going Home");
      delay(1000);
      go_home();
      //position = -19.6; //PMT Gain Correction Stand
      position = 0.0;
      Serial.println("Home (0.0 mm)");
      Serial.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    }
  }
  else if (digitalRead(RUN)) {
    Serial.println("Starting Run");
    run();
    Serial.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
  }
  int place = ((abs(position)) / 60.0) + 1;
  if (position == -19.6){
    place = 0;
  }
  //Serial.println(place);
  /*switch (place) {
    case 1:
      Serial.println("PMT #5");
      break;
    case 2:
      Serial.println("PMT #4");
      break;
    case 3:
      Serial.println("PMT #3");
      break;
    case 4:
      Serial.println("PMT #2");
      break;
    case 5:
      Serial.println("PMT #1");
    //default:
    //  Serial.println(position); 
}
  */
}

int dist_to_steps(double distance) {
  steps = round(distance/conversion);
  return steps;
}

void power_on() {
  //Serial.println("Power ON");
  //analogWrite(ENA, 255);
  analogWrite(ENB, 200);
}

void power_off() {
  //analogWrite(ENA, 0);
  analogWrite(ENB, 0);
  //Serial.println("Power OFF");
}

void move_forward(double dist) {
  steps = dist_to_steps(dist);
  //Serial.print("moving ");
  //Serial.println(steps);
  power_on();
  myStepper.step(-steps);
  power_off();
  if (digitalRead(HOME) == 1) {
    //position = -19.6; //PMT Gain Correction Stand
    position = 0.0; //S800 MPGD-DC
  }
  else {
    position = position - (steps*conversion);
  }
  Serial.print("Now at ");
  Serial.print(position);
  Serial.println(" mm.");
}

void move_backward(double dist) {
  steps = dist_to_steps(dist);
  power_on();
  myStepper.step(steps);
  power_off();
  position = position + (steps*conversion);
  Serial.print("Now at  ");
  Serial.print(position);
  Serial.println(" mm.");
}

void home() {
  Serial.println("Home Limit!");
  power_off();
  docked = 1;
  //position = -20.6;
  delay(100);
}

void end() {
  //power_off();
  position = 600.0;
  Serial.print("Reached Extreme Limit! (Position: ");
  Serial.print(position);
  Serial.println(")");
}

void run() {
  move_backward(1.0); 
  delay(2000);
  position = position - 1.0;
  move_backward(first_position);
  delay(2000);
  //for (int i =1; i < 5; i++) {
  //  move_backward(PMT_separation);
  //  delay(time_at_PMT);
  //}
}

void go_home(){
  while (docked != 1) {
      move_backward(10.0);
  }
  power_off();
  //Serial.println("Power OFF");
}
  







